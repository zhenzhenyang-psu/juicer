#!/usr/bin/awk -f
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Script to read in individual split outputs and put stats together
# Juicer version 1.5

/^[0-9]/ {
	a1+=$1; # total reads
	a2+=$4; # two_fragments
	a3+=$5; # gt2fragments
	a4+=$3; # singleton
	a5+=$2; # unmapped
}

END{
	printf("Total sequenced Reads:  %'d\nUnmapped: %'d (%0.2f%)\nSingleton: %'d (%0.2f%)\nTwo-fragment reads: %'d (%0.2f%)\ngt2fragment reads: %'d (%0.2f%)\nChimeric reads: %'d (%0.2f%)\nAlignable (Singleton+Chimeric): %'d (%0.2f%)\n", a1, a5, a5*100/a1, a4, a4*100/a1, a2, a2*100/a1, a3, a3*100/a1, a2+a3, (a2+a3)*100/a1,  a2+a3+a4, (a2+a3+a4)*100/a1);

    #printf("Total sequenced Reads:  %'d\nUnmapped: %'d (%0.2f%)\nSingleton: %'d (%0.2f%)\nTwo-fragment reads: %'d (%0.2f%)\ngt2fragment reads: %'d (%0.2f%)\nChimeric reads: %'d (%0.2f%)\nLigation Motif Present: %'d (%0.2f%)\nAlignable (Singleton+Chimeric): %'d (%0.2f%)\n", a1, a5, a5*100/a1, a4, a4*100/a1, a2, a2*100/a1, a3, a3*100/a1, a2+a3, (a2+a3)*100/a1, a6, a6*100/a1, a2+a3+a4, (a2+a3+a4)*100/a1);
}
