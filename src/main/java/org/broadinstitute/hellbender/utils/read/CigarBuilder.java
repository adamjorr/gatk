package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class CigarBuilder {

    private final List<CigarElement> cigarElements = new ArrayList<>();

    // track the last operator so we can merge consecutive elements with the same operator
    // for example, adding 3M and 4M is equivalent to adding 7M
    // also we ignore leading deletions so for example 10S + 5D = 10S
    private CigarOperator lastOperator = null;

    private Section section = Section.LEFT_HARD_CLIP;

    public CigarBuilder() { }

    public CigarBuilder add(final CigarElement element) {
        final CigarOperator operator = element.getOperator();

        if ((lastOperator == null || lastOperator.isClipping()) && !operator.consumesReadBases() && !operator.isClipping()) {
            // skip a deletion after clipping ie at the beginning of the read
            return this;
        }

        advanceSectionAndValidateCigarOrder(operator);

        // merge consecutive elements with the same operator
        if (operator == lastOperator) {
            final int n = cigarElements.size() - 1;
            cigarElements.set(n, new CigarElement(cigarElements.get(n).getLength() + element.getLength(), operator));
        } else {
            if (lastOperator == null) {
                cigarElements.add(element);
                lastOperator = operator;
            } else if (operator.isClipping() && !lastOperator.consumesReadBases() && !lastOperator.isClipping()) {
                // if we have just start clipping on the right and realize the last operator was a deletion, remove it
                cigarElements.set(cigarElements.size() - 1, element);
                lastOperator = operator;
            } else {
                cigarElements.add(element);
                lastOperator = operator;
            }
        }

        return this;
    }

    public CigarBuilder addAll(final Collection<CigarElement> elements) {
        for (final CigarElement element : elements) {
            add(element);
        }
        return this;
    }

    public Cigar make() {
        Utils.validate(section != Section.LEFT_HARD_CLIP && section != Section.LEFT_SOFT_CLIP, () -> "cigar " + cigarElements + " is completely clipped");
        return new Cigar(cigarElements);
    }

    private enum Section {LEFT_HARD_CLIP, LEFT_SOFT_CLIP, MIDDLE, RIGHT_SOFT_CLIP, RIGHT_HARD_CLIP}

    // validate that cigar structure is hard clip, soft clip, unclipped, soft clip, hard clip
    private void advanceSectionAndValidateCigarOrder(CigarOperator operator) {
        if (operator == CigarOperator.HARD_CLIP) {
            Utils.validate(section != Section.LEFT_SOFT_CLIP, "cigar is completely clipped");
            if (section == Section.MIDDLE || section == Section.RIGHT_SOFT_CLIP) {
                section = Section.RIGHT_HARD_CLIP;
            }
        } else if (operator == CigarOperator.SOFT_CLIP) {
            Utils.validate(section != Section.RIGHT_HARD_CLIP, "cigar has already reached its right hard clip");
            if (section == Section.LEFT_HARD_CLIP) {
                section = Section.LEFT_SOFT_CLIP;
            } else if(section == Section.MIDDLE) {
                section = Section.RIGHT_SOFT_CLIP;
            }
        } else {
            Utils.validate(section != Section.RIGHT_SOFT_CLIP && section != Section.RIGHT_HARD_CLIP, "cigar has already reached right clip");
            if (section == Section.LEFT_HARD_CLIP || section == Section.LEFT_SOFT_CLIP) {
                section = Section.MIDDLE;
            }
        }
    }
}