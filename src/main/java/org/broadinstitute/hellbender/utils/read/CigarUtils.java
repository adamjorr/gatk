package org.broadinstitute.hellbender.utils.read;

import com.google.common.collect.Lists;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;

public final class CigarUtils {

    // used in the bubble state machine to apply Smith-Waterman to the bubble sequence
    // these values were chosen via optimization against the NA12878 knowledge base
    public static final SWParameters NEW_SW_PARAMETERS = new SWParameters(200, -150, -260, -11);

    // In Mutect2 and HaplotypeCaller reads are realigned to their *best* haplotypes, which is very different from a generic alignment.
    // The {@code NEW_SW_PARAMETERS} penalize a substitution error more than an indel up to a length of 9 bases!
    // Suppose, for example, that a read has a single substitution error, say C -> T, on its last base.  Those parameters
    // would prefer to extend a deletion until the next T on the reference is found in order to avoid the substitution, which is absurd.
    // Since these parameters are for aligning a read to the biological sequence we believe it comes from, the parameters
    // we choose should correspond to sequencer error.  They *do not* have anything to do with the prevalence of true variation!
    public static final SWParameters ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS = new SWParameters(10, -15, -30, -5);

    private static final String SW_PAD = "NNNNNNNNNN";

    private CigarUtils(){}

    /**
     * Inverts the order of the operators in the cigar.
     * Eg 10M1D20M -> 20M1D10M
     */
    public static Cigar invertCigar (final Cigar cigar) {
        Utils.nonNull(cigar);
        final List<CigarElement>  els = new ArrayList<>(cigar.getCigarElements());
        Collections.reverse(els);
        return new Cigar(els);
    }

    /**
     * Compute the number of reference bases between the start (inclusive) and end (exclusive) cigar elements.
     * Reference bases are counted as the number of positions w.r.t. the reference spanned by an alignment to that reference.
     * For example original position = 10. cigar: 2M3I2D1M. If you remove the 2M the new starting position is 12.
     * If you remove the 2M3I it is still 12. If you remove 2M3I2D (not reasonable cigar), you will get position 14.
     */
    public static int countRefBasesAndClips(final Cigar cigar, final int cigarStartIndex, final int cigarEndIndex){
        return countRefBasesAndMaybeAlsoClips(cigar, cigarStartIndex, cigarEndIndex, true, true);
    }

    public static int countRefBasesAndSoftClips(final Cigar cigar, final int cigarStartIndex, final int cigarEndIndex){
        return countRefBasesAndMaybeAlsoClips(cigar, cigarStartIndex, cigarEndIndex, true, false);
    }

    private static int countRefBasesAndMaybeAlsoClips(final Cigar cigar, final int cigarStartIndex, final int cigarEndIndex, final boolean includeSoftClips, final boolean includeHardClips) {
        final List<CigarElement> elems = Utils.nonNull(cigar).getCigarElements();
        Utils.validateArg(cigarStartIndex >= 0 && cigarEndIndex <= elems.size() && cigarStartIndex <= cigarEndIndex, () -> "invalid index:" + 0 + " -" + elems.size());
        final Predicate<CigarOperator> consumesRefBases = CigarOperator::consumesReferenceBases;
        final Predicate<CigarOperator> clippingPredicate = op -> (includeSoftClips && op == CigarOperator.SOFT_CLIP) || (includeHardClips && op == CigarOperator.HARD_CLIP);
        final Predicate<CigarOperator> countOperator = consumesRefBases.or(clippingPredicate);

        return elems.subList(cigarStartIndex, cigarEndIndex).stream()
                .filter(elem -> countOperator.test(elem.getOperator()))
                .mapToInt(CigarElement::getLength)
                .sum();
    }

    /**
     * Removes all clipping operators from the cigar.
     */
    public static Cigar removeClipsAndPadding(final Cigar cigar) {
        Utils.nonNull(cigar, "cigar is null");
        final List<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());
        for ( final CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOperator(ce.getOperator()) ) {
                elements.add(ce);
            }
        }
        return new Cigar(elements);
    }

    private static boolean isClipOperator(final CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    private static boolean isClipOperator(final CigarElement el) {
        return isClipOperator(el.getOperator());
    }

    /**
     * Given a cigar1 and a read with cigar2,
     * this method creates cigar3 such that it has flanking clip operators from cigar2
     * and it has all operators from cigar1 in the middle.
     *
     * In other words if:
     * cigar2 = leftClip2 + noclips2 + rightClip2
     *
     * then
     * cigar3 = leftClip2 + cigar1 + rightClip2
     */
    public static Cigar reclipCigar(final Cigar cigar, final GATKRead read) {
        Utils.nonNull(cigar, "cigar");
        Utils.nonNull(read, "read");

        final List<CigarElement> elements = new ArrayList<>();
        int i = 0;
        final Cigar readCigar = read.getCigar();
        final int n = readCigar.numCigarElements();
        final List<CigarElement> readEls = readCigar.getCigarElements();

        //copy head clips
        while ( i < n && isClipOperator(readEls.get(i)) ) {
            elements.add(readEls.get(i));
            i++;
        }

        elements.addAll(cigar.getCigarElements());

        //skip over non-clips
        i++;
        while ( i < n && !isClipOperator(readEls.get(i)) ) {
            i++;
        }

        //copy tail clips
        while ( i < n && isClipOperator(readEls.get(i)) ) {
            elements.add(readEls.get(i));
            i++;
        }

        return new Cigar(elements);
    }

    /**
     * Returns whether the list has any N operators.
     */
    public static boolean containsNOperator(final List<CigarElement> cigarElements) {
        return Utils.nonNull(cigarElements).stream().anyMatch(el -> el.getOperator() == CigarOperator.N);
    }

    /**
     * A good Cigar object obeys the following rules:
     *  - is valid as per SAM spec {@link Cigar#isValid(String, long)}.
     *  - has no consecutive I/D elements
     *  - does not start or end with deletions (with or without preceding clips).
     */
    public static boolean isGood(final Cigar c) {
        Utils.nonNull(c, "cigar is null");

        //Note: the string comes from the SAMRecord so it must be a wellformed CIGAR (that is, in "\*|([0-9]+[MIDNSHPX=])+" as per SAM spec).
        //We don't have to check that
        if (c.isValid(null, -1) != null){  //if it's invalid, then it's not good
            return false;
        }
        final List<CigarElement> elems = c.getCigarElements();
        return !(hasConsecutiveIndels(elems) || startsOrEndsWithDeletionIgnoringClips(elems));
    }

    /**
     * Checks if cigar has consecutive I/D elements.
     */
    private static boolean hasConsecutiveIndels(final List<CigarElement> elems) {
        boolean prevIndel = false;
        for (final CigarElement elem : elems) {
            final CigarOperator op = elem.getOperator();
            final boolean isIndel = (op == CigarOperator.INSERTION || op == CigarOperator.DELETION);
            if (prevIndel && isIndel) {
                return true;
            }
            prevIndel = isIndel;
        }
        return false;
    }

    /**
     * Checks if cigar starts with a deletion (ignoring any clips at the beginning).
     */
    private static boolean startsOrEndsWithDeletionIgnoringClips(final List<CigarElement> elems) {
        CigarOperator lastOp = null;
        for (final CigarElement elem : elems) {
               final CigarOperator op = elem.getOperator();
               if ((lastOp != null && lastOp.isClipping() && op == CigarOperator.DELETION)
                   || (lastOp == CigarOperator.DELETION && op.isClipping())) {
                   return true;
               }
               lastOp = op;
        }
        return false;
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param aligner
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public static Cigar calculateCigar(final byte[] refSeq, final byte[] altSeq, final SmithWatermanAligner aligner, final SWOverhangStrategy strategy) {
        Utils.nonNull(refSeq, "refSeq");
        Utils.nonNull(altSeq, "altSeq");
        if ( altSeq.length == 0 ) {
            // horrible edge case from the unit tests, where this path has no bases
            return new Cigar(Collections.singletonList(new CigarElement(refSeq.length, CigarOperator.D)));
        }

        //Note: this is a performance optimization.
        // If two strings are equal (a O(n) check) then it's trivial to get CIGAR for them.
        // Furthermore, if their lengths are equal and their element-by-element comparison yields two or fewer mismatches
        // it's also a trivial M-only CIGAR, because in order to have equal length one would need at least one insertion and
        // one deletion, in which case two substitutions is a better alignment.
        if (altSeq.length == refSeq.length){
            int mismatchCount = 0;
            for (int n = 0; n < refSeq.length && mismatchCount <= 2; n++) {
                mismatchCount += (altSeq[n] == refSeq[n] ? 0 : 1);
            }
            if (mismatchCount <= 2) {
                final Cigar matching = new Cigar();
                matching.add(new CigarElement(refSeq.length, CigarOperator.MATCH_OR_MISMATCH));
                return matching;
            }
        }

        final Cigar nonStandard;

        final String paddedRef = SW_PAD + new String(refSeq) + SW_PAD;
        final String paddedPath = SW_PAD + new String(altSeq) + SW_PAD;
        final SmithWatermanAlignment alignment = aligner.align(paddedRef.getBytes(), paddedPath.getBytes(), NEW_SW_PARAMETERS, strategy);

        if ( isSWFailure(alignment) ) {
            return null;
        }


        // cut off the padding bases
        final int baseStart = SW_PAD.length();
        final int baseEnd = paddedPath.length() - SW_PAD.length() - 1; // -1 because it's inclusive
        nonStandard = AlignmentUtils.trimCigarByBases(alignment.getCigar(), baseStart, baseEnd);

        if ( nonStandard.getReferenceLength() != refSeq.length ) {
            nonStandard.add(new CigarElement(refSeq.length - nonStandard.getReferenceLength(), CigarOperator.D));
        }

        // finally, return the cigar with all indels left aligned
        return leftAlignCigarSequentially(nonStandard, refSeq, altSeq, 0, 0);
    }

    /**
     * Make sure that the SW didn't fail in some terrible way, and throw exception if it did
     */
    private static boolean isSWFailure(final SmithWatermanAlignment alignment) {
        // check that the alignment starts at the first base, which it should given the padding
        if ( alignment.getAlignmentOffset() > 0 ) {
            return true;
        }

        // check that we aren't getting any S operators (which would be very bad downstream)
        for ( final CigarElement ce : alignment.getCigar().getCigarElements() ) {
            if ( ce.getOperator() == CigarOperator.S )
                return true;
            // soft clips at the end of the alignment are really insertions
        }

        return false;
    }

    /**
     * Left align the given cigar sequentially. This is needed because AlignmentUtils doesn't accept cigars with more than one indel in them.
     * This is a target of future work to incorporate and generalize into AlignmentUtils for use by others.
     * @param cigar     the cigar to left align
     * @param refSeq    the reference byte array
     * @param readSeq   the read byte array
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @return          the left-aligned cigar
     */
    public static Cigar leftAlignCigarSequentially(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        Utils.nonNull(cigar, "cigar null");
        Utils.nonNull(refSeq, "refSeq null");
        Utils.nonNull(readSeq, "readSeq null");

        final Cigar cigarToReturn = new Cigar();
        Cigar cigarToAlign = new Cigar();
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            final CigarElement ce = cigar.getCigarElement(i);
            if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
                cigarToAlign.add(ce);
                final Cigar leftAligned = AlignmentUtils.leftAlignSingleIndel(cigarToAlign, refSeq, readSeq, refIndex, readIndex, false);
                for ( final CigarElement toAdd : leftAligned.getCigarElements() ) { cigarToReturn.add(toAdd); }
                refIndex += cigarToAlign.getReferenceLength();
                readIndex += cigarToAlign.getReadLength();
                cigarToAlign = new Cigar();
            } else {
                cigarToAlign.add(ce);
            }
        }
        if( !cigarToAlign.isEmpty() ) {
            for( final CigarElement toAdd : cigarToAlign.getCigarElements() ) {
                cigarToReturn.add(toAdd);
            }
        }

        final Cigar result = AlignmentUtils.consolidateCigar(cigarToReturn);
        Utils.validate(result.getReferenceLength() == cigar.getReferenceLength(),
                () -> "leftAlignCigarSequentially failed to produce a valid CIGAR.  Reference lengths differ.  Initial cigar " + cigar + " left aligned into " + result);
        return result;
    }

    /**
     * Returns the length of the original read considering all clippings based on this cigar.
     * <p>
     *     The result of applying this method on a empty cigar is zero.
     * </p>
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if the input {@code cigar} is {@code null}.
     * @return 0 or greater.
     */
    public static int countUnclippedReadBases(final Cigar cigar) {
        Utils.nonNull(cigar, "the input cigar cannot be null");
        return cigar.getCigarElements().stream()
                .filter(ce -> {
                    final CigarOperator operator = ce.getOperator();
                    return operator.isClipping() || operator.consumesReadBases();
                })
                .mapToInt(CigarElement::getLength)
                .sum();
    }

    private static int countClippedBases(final Cigar cigar, final ClippingTail tail, final boolean includeSoftClips, final boolean includeHardClips) {
        Utils.nonNull(cigar);
        if (cigar.numCigarElements() < 2) {
            return 0;
        } else {
            Utils.validate(includeHardClips || includeSoftClips, "no clips chosen");
            final Predicate<CigarOperator> pred = !includeHardClips ? op -> op == CigarOperator.S :
                    (includeSoftClips ? op -> op.isClipping() : op -> op == CigarOperator.H);
            int result = 0;
            for (final CigarElement elem : tail == ClippingTail.LEFT_TAIL ? cigar : Lists.reverse(cigar.getCigarElements())) {
                if (!elem.getOperator().isClipping()) {
                    return result;
                } else if (pred.test(elem.getOperator())) {
                    result += elem.getLength();
                }
            }
        }

        throw new IllegalArgumentException("the input cigar only contains clips!");
    }

    /**
     * Total number of bases clipped on the left/head side of the cigar.
     *
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if {@code cigar} is {@code null}.
     * @return 0 or greater.
     */
    public static int countLeftClippedBases(final Cigar cigar) {
        return countClippedBases(cigar, ClippingTail.LEFT_TAIL, true, true);
    }

    /**
     * Returns the number of based hard-clipped to the left/head of the cigar.
     *
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if {@code cigar} is {@code null}.
     * @return 0 or greater.
     */
    public static int countLeftHardClippedBases(final Cigar cigar) {
        return countClippedBases(cigar, ClippingTail.LEFT_TAIL, false, true);
    }

    /**
     * Returns the number of based hard-clipped to the right/tail of the cigar.
     *
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if {@code cigar} is {@code null}.
     * @return 0 or greater.
     */
    public static int countRightHardClippedBases(final Cigar cigar) {
        return countClippedBases(cigar, ClippingTail.RIGHT_TAIL, false, true);
    }

    /**
     * Total number of bases clipped (soft or hard) on the right/tail side of the cigar.
     *
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if {@code cigar} is {@code null}
     * @return 0 or greater.
     */
    public static int countRightClippedBases(final Cigar cigar) {
        return countClippedBases(cigar, ClippingTail.RIGHT_TAIL, true, true);
    }

    public static int countAlignedBases(final Cigar cigar ) {
        return Utils.nonNull(cigar).getCigarElements().stream()
                .filter(cigarElement -> cigarElement.getOperator().isAlignment())
                .mapToInt(CigarElement::getLength)
                .sum();
    }

    /**
     * replace soft clips (S) with match (M) operators, merging any consecutive M's that result eg 10S10M -> 20M.
     */
    public static Cigar revertSoftClips(final Cigar originalCigar) {
        final CigarBuilder builder = new CigarBuilder();
        for (final CigarElement element : originalCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.SOFT_CLIP) {
                builder.add(new CigarElement(element.getLength(), CigarOperator.MATCH_OR_MISMATCH));
            } else {
                builder.add(element);
            }
        }

        return builder.make();
    }

    /**
     * Given a cigar string, soft clip up to leftClipEnd and soft clip starting at rightClipBegin
     * @param start initial index to clip within read bases, inclusive
     * @param stop final index to clip within read bases exclusive
     * @param clippingOperator      type of clipping -- must be either hard clip or soft clip
     */
    public static Cigar clipCigar(final Cigar cigar, final int start, final int stop, CigarOperator clippingOperator) {
        Utils.validateArg(clippingOperator.isClipping(), "Not a clipping operator");
        final boolean clipLeft = start == 0;

        final CigarBuilder newCigar = new CigarBuilder();

        int elementStart = 0;
        for (final CigarElement element : cigar.getCigarElements()) {
            final CigarOperator operator = element.getOperator();
            // copy hard clips
            if (operator == CigarOperator.HARD_CLIP) {
                newCigar.add(new CigarElement(element.getLength(), element.getOperator()));
                continue;
            }
            final int elementEnd = elementStart + (operator.consumesReadBases() ? element.getLength() : 0);

            // element precedes start or follows end of clip, copy it to new cigar
            if (elementEnd <= start || elementStart >= stop) {
                newCigar.add(new CigarElement(element.getLength(), operator));
            } else {    // otherwise, some or all of the element is soft-clipped
                final int unclippedLength = clipLeft ? elementEnd - stop : start - elementStart;
                final int clippedLength = element.getLength() - unclippedLength;

                if (unclippedLength <= 0) { // totally clipped
                    if (operator.consumesReadBases()) {
                        newCigar.add(new CigarElement(element.getLength(), clippingOperator));
                    }
                } else if (clipLeft) {
                    newCigar.add(new CigarElement(clippedLength, clippingOperator));
                    newCigar.add(new CigarElement(unclippedLength, operator));
                } else {
                    newCigar.add(new CigarElement(unclippedLength, operator));
                    newCigar.add(new CigarElement(clippedLength, clippingOperator));
                }
            }
            elementStart = elementEnd;
        }

        return newCigar.make();
    }

    /**
     * How many bases to the right does a read's alignment start shift given its cigar and the number of left soft clips
     */
    public static int alignmentStartShift(final Cigar cigar, final int numClipped) {
        int refBasesClipped = 0;

        int elementStart = 0;   // this and elementEnd are indices in the read's bases
        for (final CigarElement element : cigar.getCigarElements()) {
            final CigarOperator operator = element.getOperator();
            // hard clips consume neither read bases nor reference bases and are irrelevant
            if (operator == CigarOperator.HARD_CLIP) {
                continue;
            }
            final int elementEnd = elementStart + (operator.consumesReadBases() ? element.getLength() : 0);

            if (elementEnd <= numClipped) {  // totally within clipped span -- this includes deletions immediately following clipping
                refBasesClipped += operator.consumesReferenceBases() ? element.getLength() : 0;
            } else if (elementStart < numClipped) { // clip in middle of element, which means the element necessarily consumes read bases
                final int clippedLength = numClipped - elementStart;
                refBasesClipped += operator.consumesReferenceBases() ? clippedLength : 0;
                break;
            }
            elementStart = elementEnd;
        }
        return refBasesClipped;
    }
}
