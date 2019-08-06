Refer to MacWilliams and Sloane "The theory of Error-Correcting Codes" for all the terms used here.

This repository contains programs for solving 3 tasks from the coding theory:
1. It builds an error-correcting code, which is capable of correcting arbitrary number of errors.
2. It codes an input message according to the previously built code and then it distorts the codeword.
3. It decodes a distorted codeword using maximum likelihood decoding.

I. For building a code it uses the Gilbert-Varshamov bound to determine an accepteable set of parameters: 
- number of check symbols,
- a codeword length,
- number of errors to correct.

Then it generates 3 files: 
- a generator matrix file, 
- a parity check matrix file,
- a standard array file.

The 1st file is used in II, the 2nd and the 3rd are used in III.

II. The generator matrix is used to code a message. 
Then an arbitrary or a particular error vector is added to the codeword.

III. The parity check matrix and the standard array are used to decode a distorted message. 