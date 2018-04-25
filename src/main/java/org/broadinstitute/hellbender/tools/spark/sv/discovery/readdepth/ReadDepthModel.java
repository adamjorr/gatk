package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.solver.SimulatedAnnealingSolver;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

public final class ReadDepthModel implements Serializable {

    public static final long serialVersionUID = 1L;
    private final Map<SimpleSVType.TYPES,List<ReadDepthCluster>> clusteredEvents;
    private ReadDepthModelParameters parameters;
    private final Logger logger = LogManager.getLogger(this.getClass());
    private long seed;

    public ReadDepthModel(final SVIntervalTree<LargeSimpleSV> eventsTree, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        this.parameters = new ReadDepthModelParameters();
        setSamplerSeed(0);
        this.clusteredEvents = clusterEvents(eventsTree, copyRatioSegmentOverlapDetector, dictionary);
    }

    private static double unscaledNormal(final double x, final double mu, final double sigma) {
        final double diff = (x - mu) / sigma;
        return Math.exp(-0.5 * diff * diff) / sigma;
    }

    private static double unscaledLogNormal(final double x, final double mu, final double sigma) {
        final double diff = (x - mu) / sigma;
        return -0.5 * diff * diff - Math.log(sigma);
    }

    public void setParameters(final ReadDepthModelParameters parameters) {
        this.parameters = parameters;
    }

    public void setSamplerSeed(final long seed) {
        this.seed = seed;
    }

    private Map<SimpleSVType.TYPES,List<ReadDepthCluster>> clusterEvents(final SVIntervalTree<LargeSimpleSV> eventsTree, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        final Set<LargeSimpleSV> visited = new HashSet<>(SVUtils.hashMapCapacity(eventsTree.size()));
        final Map<SimpleSVType.TYPES,List<ReadDepthCluster>> clusteredEvents = new HashMap<>(SVUtils.hashMapCapacity(eventsTree.size()));
        final Iterator<SVIntervalTree.Entry<LargeSimpleSV>> iter = eventsTree.iterator();
        final Set<LargeSimpleSV> cluster = new HashSet<>(SVUtils.hashMapCapacity(eventsTree.size()));
        while (iter.hasNext()) {
            final SVIntervalTree.Entry<LargeSimpleSV> entry = iter.next();
            final LargeSimpleSV event = entry.getValue();
            if (visited.contains(event)) continue;
            cluster.clear();
            cluster.add(event);
            int oldClusterSize;
            do {
                final List<LargeSimpleSV> oldCluster = new ArrayList<>(cluster); //Have to copy the set to avoid ConcurrentModificationException
                oldClusterSize = cluster.size();
                for (final LargeSimpleSV clusterEvent : oldCluster) {
                    cluster.addAll(Utils.stream(eventsTree.overlappers(clusterEvent.getInterval()))
                            .map(SVIntervalTree.Entry::getValue)
                            .filter(overlapper -> overlapper.getType() == event.getType()) //Cluster by event type
                            .collect(Collectors.toList()));
                }
            } while (cluster.size() > oldClusterSize);
            visited.addAll(cluster);

            final List<ReadDepthEvent> newCluster = new ArrayList<>(cluster.size());
            int eventId = 0;
            final Iterator<LargeSimpleSV> newClusterIter = cluster.iterator();
            while (newClusterIter.hasNext()) {
                newCluster.add(new ReadDepthEvent(eventId, newClusterIter.next()));
                eventId++;
            }
            clusteredEvents.putIfAbsent(event.getType(), new ArrayList<>());
            clusteredEvents.get(event.getType()).add(new ReadDepthCluster(newCluster, copyRatioSegmentOverlapDetector, dictionary));
        }
        return clusteredEvents;
    }

    public Collection<ReadDepthEvent> getEvents() {
        return clusteredEvents.values().stream().flatMap(List::stream).flatMap(entry -> entry.getEventsList().stream()).collect(Collectors.toList());
    }

    private static double metropolisHastingsTemperatureSchedule(final int step) {
        return 1.0;
    }

    public Tuple2<Double,List<ReadDepthEvent>> solve(final JavaSparkContext ctx) {
        final JavaRDD<ReadDepthCluster> clusterRdd = ctx.parallelize(clusteredEvents.values().stream().flatMap(List::stream).collect(Collectors.toList()));
        final List<Tuple2<ReadDepthCluster,Double>> result = solveStateMaximumPosterior(clusterRdd, parameters, seed);
        //solveParameterMaximumPosterior(result.stream().map(Tuple2::_1).collect(Collectors.toList()), parameters, seed, logger);
        final double energy = result.stream().mapToDouble(Tuple2::_2).sum();
        final List<ReadDepthEvent> events = result.stream().map(Tuple2::_1).map(ReadDepthCluster::getEventsList).flatMap(List::stream).collect(Collectors.toList());
        logger.info("Done: f = " + energy);
        return new Tuple2<>(energy,events);
    }

    public static double solveParameterMaximumPosterior(final List<ReadDepthCluster> clusters, final ReadDepthModelParameters parameters, final long seed, final Logger logger) {
        final int size = parameters.getParameters().length;
        final NormalDistribution standardNormal = new NormalDistribution(0, 1);
        standardNormal.reseedRandomGenerator(seed);
        final Function<double[],Double> energyFunction = x -> clusters.stream().mapToDouble(cluster -> rejectDupAcceptDelEnergy(x, cluster, seed)).sum();
        final Supplier<double[]> sampler = () -> parameterStepSampler(standardNormal);
        final double[] lowerBound = ReadDepthModelParameters.getDefaultLowerBounds();
        final double[] upperBound = ReadDepthModelParameters.getDefaultUpperBounds();
        final SimulatedAnnealingSolver parameterSolver = new SimulatedAnnealingSolver(size, energyFunction, sampler, lowerBound, upperBound);

        final double[] x0 = parameters.getParameters();
        final double finalEnergy = parameterSolver.solve(x0, 3, 0);
        parameters.setParameters(parameterSolver.getSolution());
        logger.info(parameters.toString());
        return finalEnergy;
    }

    private static double rejectDupAcceptDelEnergy(final double[] x, final ReadDepthCluster cluster, final long seed) {
        double err = 0;
        final ReadDepthModelParameters parameters = new ReadDepthModelParameters();
        parameters.setParameters(x);
        solveStateMaximumPosterior(cluster, parameters, seed);
        for (final ReadDepthEvent event : cluster.getEventsList()) {
            final int id = event.getId();
            final SimpleSVType.TYPES type = cluster.getType();
            if (type == SimpleSVType.TYPES.DEL) {
                err += 1 - Math.min(1, x[id]);
            } else if (type == SimpleSVType.TYPES.DUP_TAND) {
                err += x[id];
            }
        }
        return err;
    }

    private static List<Tuple2<ReadDepthCluster,Double>> solveStateMaximumPosterior(final JavaRDD<ReadDepthCluster> clusters, final ReadDepthModelParameters parameters, final long seed) {
        return clusters.map(cluster -> solveStateMaximumPosterior(cluster, parameters, seed)).collect();
    }

    private static double[] stateStepSampler(final int size, final Random random) {
        /*final double[] sample = standardNormal.sample(size);
        for (int i = 0; i < sample.length; i++) {
            sample[i] *= 0.1;
        }
        return sample;*/
        final double[] sample = new double[size];
        boolean isZero = true;
        while (isZero) {
            for (int i = 0; i < size; i++) {
                sample[i] = random.nextInt(3) - 1;
                if (sample[i] != 0) {
                    isZero = false;
                }
            }
        }
        return sample;
    }

    public static Tuple2<ReadDepthCluster,Double> solveStateMaximumPosterior(final ReadDepthCluster cluster, final ReadDepthModelParameters parameters, final long seed) {
        if (cluster.getEventsList().size() == 1) {
            return solveStateMaximumPosteriorBF(cluster, parameters);
        }
        return solveStateMaximumPosteriorSA(cluster, parameters, seed);
    }

    private static Tuple2<ReadDepthCluster,Double> solveStateMaximumPosteriorBF(final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        final int numEvents = cluster.getEventsTree().size();
        if (numEvents != 1) {
            throw new GATKException("Brute force only supports one event");
        }
        final double[] x = new double[numEvents];
        double minX = 0;
        double minEnergy = -computeLogPosterior(x, cluster, parameters);
        final int maxPloidy = (int) parameters.getParameter(ReadDepthModelParameters.ParameterEnum.MAX_PLOIDY);
        for (int i = 1; i <= maxPloidy; i++) {
            x[0] = i;
            final double energy = -computeLogPosterior(x, cluster, parameters);
            if (i == 0 || energy < minEnergy) {
                minX = i;
                minEnergy = energy;
            }
        }
        cluster.getEventsList().get(0).setState(minX);
        return new Tuple2<>(cluster, minEnergy);
    }

    private static Tuple2<ReadDepthCluster,Double> solveStateMaximumPosteriorSA(final ReadDepthCluster cluster, final ReadDepthModelParameters parameters, final long seed) {
        final List<ReadDepthEvent> events = cluster.getEventsList();
        final Function<double[], Double> energyFunction = x -> -computeLogPosterior(x, cluster, parameters);
        final int size = events.size();
        final Random random = new Random(seed);
        final Supplier<double[]> sampler = () -> stateStepSampler(size, random);
        final double[] lowerBound = new double[size];
        final double[] upperBound = new double[size];
        Arrays.fill(lowerBound, 0);
        final int maxPloidy = (int) parameters.getParameter(ReadDepthModelParameters.ParameterEnum.MAX_PLOIDY);
        Arrays.fill(upperBound, maxPloidy);
        final SimulatedAnnealingSolver sa = new SimulatedAnnealingSolver(size, energyFunction, sampler, lowerBound, upperBound);

        final double[] x0 = new double[size];
        for (int i = 0; i < events.size(); i++) {
            x0[i] = events.get(i).getState();
        }
        final int numSteps = (int) Math.min(100000, Math.pow(10, size));
        //final double finalEnergy = sa.solve(x0, numSteps, 1000);
        final double finalEnergy = sa.solve(x0, numSteps, ReadDepthModel::metropolisHastingsTemperatureSchedule, 1000);

        final double[] minEnergyState = sa.getSolution();
        for (int i = 0; i < events.size(); i++) {
            events.get(i).setState(minEnergyState[i]);
        }
        return new Tuple2<>(cluster, finalEnergy);
    }

    private static double computeLogPosterior(final double[] x, final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        return computeLogLikelihood(x, cluster, parameters) + computeLogPrior(x, cluster, parameters);
    }

    private static double computeLogLikelihood(final double[] x, final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        return cluster.getCopyNumberInfo().stream().mapToDouble(entry -> computeCopyNumberLikelihood(x, entry, parameters)).sum()
                + cluster.getEventsList().stream().mapToDouble(event -> computeReadEvidenceLikelihood(x, event, parameters)).sum()
                + computeCallDistanceLikelihood(x, cluster, parameters);
    }

    private static double computeLogPrior(final double[] x, final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        return cluster.getCopyNumberInfo().stream().mapToDouble(entry -> computePloidyPrior(x, entry, cluster.getType(), parameters)).sum();
    }

    private static double computeCallDistanceLikelihood(final double[] x, final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        //final double meanInsertSize = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.MEAN_INSERT_SIZE);
        //final double callDistancePseudocount = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_DISTANCE_PSEUDOCOUNT);
        double total = 0;
        final List<ReadDepthEvent> events = cluster.getEventsList();
        final List<Tuple2<Double, Double>> nearestCallDistances = cluster.getNearestCallDistances();
        for (int i = 0; i < events.size(); i++) {
            final double mean = x[i] == 0 ? 6 : 4;
            final double std = x[i] == 0 ? 1000: 1000; // 1 : 1; //TODO
            total += events.get(i).getEvent().getSize() * (unscaledLogNormal(nearestCallDistances.get(i)._1, mean, std) + unscaledLogNormal(nearestCallDistances.get(i)._2, mean, std));
            //final double std = meanInsertSize / Math.max(callDistancePseudocount, callDistancePseudocount + Math.min(x[i], 1));
            //total += events.get(i).getEvent().getSize() * (unscaledLogNormal(nearestCallDistances.get(i)._1, 0, std) + unscaledLogNormal(nearestCallDistances.get(i)._2, 0, std));
        }
        return total;
    }

    private static double computeReadEvidenceLikelihood(final double[] x, final ReadDepthEvent event, final ReadDepthModelParameters parameters) {
        final double singleCopyDepth = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.COPY_NEUTRAL_DEPTH) * 0.5;
        final double expectedReadEvidenceFraction = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.EXPECTED_READ_EVIDENCE_FRACTION);
        final double expectedReadEvidenceStd = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.EXPECTED_READ_EVIDENCE_STD);
        final double expectedEvidence = x[event.getId()] * singleCopyDepth * expectedReadEvidenceFraction;
        final double sigma = expectedReadEvidenceStd * singleCopyDepth;
        return event.getEvent().getSize() * (unscaledLogNormal(event.getEvent().getReadPairEvidence(), expectedEvidence, sigma)
                + unscaledLogNormal(event.getEvent().getSplitReadEvidence(), expectedEvidence, sigma));
    }

    private static double computeCopyNumberLikelihood(final double[] x, final Tuple2<List<OverlapInfo>, Double> entry, final ReadDepthModelParameters parameters) {
        final List<OverlapInfo> overlapInfoList = entry._1;
        final double calledCopyNumber = entry._2;
        final double sigma = calledCopyNumber == 0 ? 10 : parameters.getParameter(ReadDepthModelParameters.ParameterEnum.COPY_NUMBER_STD); //TODO different distribution for CN 0?
        return overlapInfoList.stream()
                .mapToDouble(info -> {
                    final double estimatedCopyNumber = 2 + info.idsAndCoefficients.stream().mapToDouble(tuple -> tuple._2 * x[tuple._1]).sum();
                    return unscaledLogNormal(estimatedCopyNumber, calledCopyNumber, sigma) * info.size;
                }).sum();
    }

    private static double computePloidyPrior(final double[] x, final Tuple2<List<OverlapInfo>, Double> entry,
                                      final SimpleSVType.TYPES type, final ReadDepthModelParameters parameters) {
        final List<OverlapInfo> overlapInfoList = entry._1;
        final double maxPloidy = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.MAX_PLOIDY);
        final double parameterConstraintStd = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.PARAMETER_CONSTRAINT_STD);
        double total = 0;
        for (final OverlapInfo info : overlapInfoList) {
            double totalPloidy = info.idsAndCoefficients.stream().mapToDouble(pair -> getPloidy(x[pair._1], type)).sum();
            final double excessPloidy;
            if (totalPloidy > maxPloidy) {
                excessPloidy = totalPloidy - maxPloidy;
            } else {
                excessPloidy = 0;
            }
            total += unscaledLogNormal(excessPloidy, 0, parameterConstraintStd) * info.size;
        }
        return total;
    }

    private static double getPloidy(final double val, final SimpleSVType.TYPES type) {
        return type == SimpleSVType.TYPES.DUP_TAND ? Math.min(val, 1) : val;
    }

    private final static class ReadDepthModelParameters implements Serializable {

        public static final long serialVersionUID = 1L;
        public enum ParameterEnum {
            SA_STEPS,
            SA_T0,
            MAX_PLOIDY,
            COPY_NEUTRAL_DEPTH,
            MEAN_INSERT_SIZE,
            PARAMETER_CONSTRAINT_STD,
            EXPECTED_READ_EVIDENCE_FRACTION,
            EXPECTED_READ_EVIDENCE_STD,
            COPY_NUMBER_STD,
            CALL_DISTANCE_PSEUDOCOUNT
        }

        public static final int DEFAULT_SA_STEPS = 10000;
        public static final double DEFAULT_SA_T0 = 10000;
        public static final double DEFAULT_MAX_PLOIDY = 2;
        public static final double DEFAULT_COPY_NEUTRAL_DEPTH = 30;
        public static final double DEFAULT_MEAN_INSERT_SIZE = 500;
        public static final double DEFAULT_PARAMETER_CONSTRAINT_STD = 0.001;
        public static final double DEFAULT_EXPECTED_READ_EVIDENCE_FRACTION = 1.0;
        public static final double DEFAULT_EXPECTED_READ_EVIDENCE_STD = 1000; //0.2;
        public static final double DEFAULT_COPY_NUMBER_STD = 1000; //0.5;
        public static final double DEFAULT_CALL_DISTANCE_PSEUDOCOUNT = 0.01;

        private double[] parameters;

        public ReadDepthModelParameters() {
            parameters = new double[ParameterEnum.values().length];
            setParameter(ParameterEnum.SA_STEPS, (double) DEFAULT_SA_STEPS);
            setParameter(ParameterEnum.SA_T0, DEFAULT_SA_T0);
            setParameter(ParameterEnum.MAX_PLOIDY, DEFAULT_MAX_PLOIDY);
            setParameter(ParameterEnum.COPY_NEUTRAL_DEPTH, DEFAULT_COPY_NEUTRAL_DEPTH);
            setParameter(ParameterEnum.MEAN_INSERT_SIZE, DEFAULT_MEAN_INSERT_SIZE);
            setParameter(ParameterEnum.PARAMETER_CONSTRAINT_STD, DEFAULT_PARAMETER_CONSTRAINT_STD);
            setParameter(ParameterEnum.EXPECTED_READ_EVIDENCE_FRACTION, DEFAULT_EXPECTED_READ_EVIDENCE_FRACTION);
            setParameter(ParameterEnum.EXPECTED_READ_EVIDENCE_STD, DEFAULT_EXPECTED_READ_EVIDENCE_STD);
            setParameter(ParameterEnum.COPY_NUMBER_STD, DEFAULT_COPY_NUMBER_STD);
            setParameter(ParameterEnum.CALL_DISTANCE_PSEUDOCOUNT, DEFAULT_CALL_DISTANCE_PSEUDOCOUNT);
        }

        public double[] getParameters() {
            return parameters;
        }

        public void setParameter(final ParameterEnum parameter, final double value) {
            parameters[parameter.ordinal()] = value;
        }

        public double getParameter(final ParameterEnum parameter) {
            return parameters[parameter.ordinal()];
        }

        public void setParameters(final double[] values) {
            Utils.nonNull(values, "Parameter vector cannot be null");
            Utils.validateArg(values.length == ParameterEnum.values().length, "Invalid parameter vector size");
            parameters = values;
        }

        public static double[] getDefaultUpperBounds() {
            final double[] params = new double[ParameterEnum.values().length];
            params[ParameterEnum.SA_STEPS.ordinal()] = DEFAULT_SA_STEPS;
            params[ParameterEnum.SA_T0.ordinal()] = DEFAULT_SA_T0;
            params[ParameterEnum.MAX_PLOIDY.ordinal()] = DEFAULT_MAX_PLOIDY;
            params[ParameterEnum.COPY_NEUTRAL_DEPTH.ordinal()] = DEFAULT_COPY_NEUTRAL_DEPTH;
            params[ParameterEnum.MEAN_INSERT_SIZE.ordinal()] = 1000;
            params[ParameterEnum.PARAMETER_CONSTRAINT_STD.ordinal()] = 0.1;
            params[ParameterEnum.EXPECTED_READ_EVIDENCE_FRACTION.ordinal()] = 2;
            params[ParameterEnum.EXPECTED_READ_EVIDENCE_STD.ordinal()] = 1;
            params[ParameterEnum.COPY_NUMBER_STD.ordinal()] = 10;
            params[ParameterEnum.CALL_DISTANCE_PSEUDOCOUNT.ordinal()] = 10;
            return params;
        }

        public static double[] getDefaultLowerBounds() {
            final double[] params = new double[ParameterEnum.values().length];
            params[ParameterEnum.SA_STEPS.ordinal()] = DEFAULT_SA_STEPS;
            params[ParameterEnum.SA_T0.ordinal()] = DEFAULT_SA_T0;
            params[ParameterEnum.MAX_PLOIDY.ordinal()] = DEFAULT_MAX_PLOIDY;
            params[ParameterEnum.COPY_NEUTRAL_DEPTH.ordinal()] = DEFAULT_COPY_NEUTRAL_DEPTH;
            params[ParameterEnum.MEAN_INSERT_SIZE.ordinal()] = 0;
            params[ParameterEnum.PARAMETER_CONSTRAINT_STD.ordinal()] = 1e-3;
            params[ParameterEnum.EXPECTED_READ_EVIDENCE_FRACTION.ordinal()] = 0;
            params[ParameterEnum.EXPECTED_READ_EVIDENCE_STD.ordinal()] = 1e-3;
            params[ParameterEnum.COPY_NUMBER_STD.ordinal()] = 1e-3;
            params[ParameterEnum.CALL_DISTANCE_PSEUDOCOUNT.ordinal()] = 1e-10;
            return params;
        }

        @Override
        public String toString() {
            final StringBuilder stringBuilder = new StringBuilder();
            for (final ParameterEnum type : ParameterEnum.values()) {
                stringBuilder.append("{");
                stringBuilder.append(type.toString());
                stringBuilder.append(": ");
                stringBuilder.append(getParameter(type));
                stringBuilder.append("}");
            }
            return stringBuilder.toString();
        }
    }

    final static double[] parameterStepSampler(final NormalDistribution standardNormal) {
        final double[] sample = standardNormal.sample(ReadDepthModelParameters.ParameterEnum.values().length);
        sample[ReadDepthModelParameters.ParameterEnum.SA_STEPS.ordinal()] = 0;
        sample[ReadDepthModelParameters.ParameterEnum.SA_T0.ordinal()] = 0;
        sample[ReadDepthModelParameters.ParameterEnum.MAX_PLOIDY.ordinal()] = 0;
        sample[ReadDepthModelParameters.ParameterEnum.COPY_NEUTRAL_DEPTH.ordinal()] = 0;
        sample[ReadDepthModelParameters.ParameterEnum.MEAN_INSERT_SIZE.ordinal()] *= 50;
        sample[ReadDepthModelParameters.ParameterEnum.PARAMETER_CONSTRAINT_STD.ordinal()] *= 0.005;
        sample[ReadDepthModelParameters.ParameterEnum.EXPECTED_READ_EVIDENCE_FRACTION.ordinal()] *= 0.1;
        sample[ReadDepthModelParameters.ParameterEnum.EXPECTED_READ_EVIDENCE_STD.ordinal()] *= 0.1;
        sample[ReadDepthModelParameters.ParameterEnum.COPY_NUMBER_STD.ordinal()] *= 0.05;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_DISTANCE_PSEUDOCOUNT.ordinal()] *= 0.01;
        return sample;
    }

    static final class OverlapInfo implements Serializable {
        public static final long serialVersionUID = 1L;
        public List<Tuple2<Integer, Integer>> idsAndCoefficients;
        public double size;

        public OverlapInfo(List<Tuple2<Integer, Integer>> idsAndCoefficients, double size) {
            this.idsAndCoefficients = idsAndCoefficients;
            this.size = size;
        }
    }

    /*
    private double probabilityGivenZygosityPrior(final double[] x, final ReadDepthCluster cluster) {
        double total = 0;
        for (int i = 0; i < r.length; i++) {
            total += unscaledLogNormal(Math.min(1, q[i]) - r[i], 0, parameters.pqDifferenceStd) * events.get(i).getSize();
        }
        return total;
    }
    */

    /*
    private double quantizationPrior(final double[] x, final int numStates, final ReadDepthCluster cluster) {
        double total = 0;
        for (int i = 0; i < x.length; i++) {
            final double size = events.get(i).getSize();
            for (int j = 0; j < numStates; j++) {
                total += unscaledLogNormal(x[i], j, 1.0) * size;
            }
        }
        return total;
    }
    */

    /*
    private double parameterBoundsPrior(final double[] x, final double xMax, final ReadDepthCluster cluster) {
        double total = 0;
        for (int i = 0; i < x.length; i++) {
            final double y;
            if (x[i] > xMax) {
                y = xMax - x[i];
            } else if (x[i] < 0) {
                y = x[i];
            } else {
                y = 0;
            }
            total += unscaledLogNormal(y, 0, parameterConstraintStd) * events.get(i).getSize();
        }
        return total;
    }
    */
}
