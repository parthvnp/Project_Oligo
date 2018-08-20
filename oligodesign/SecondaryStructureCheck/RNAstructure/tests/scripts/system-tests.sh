verifyInputFiles $SINGLESEQ && echo "GOOD" || echo "BAD"

verifyInputFiles banana.txt && echo "BAD" || echo "GOOD"

(verifyInputFiles banana.txt -r; echo "BAD") && echo "BAD" || echo "GOOD"
