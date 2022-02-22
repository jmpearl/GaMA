# runtests.sh
LOGFILE=log.txt

matlab -nodesktop -nosplash -minimize -wait -logfile "$LOGFILE" -r 'scripts/Tests/Test_All';
CODE=$?

cat "$LOGFILE"

exit $CODE
