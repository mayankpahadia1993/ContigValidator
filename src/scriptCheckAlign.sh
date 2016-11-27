# errorRates=( 0 0.01 0.015 0.020 0.025 0.030 0.035)
errorRates=( 0.02 0.025 0.03 0.035)
# errorRates=( 0.02)
for i in "${errorRates[@]}"
do
	echo $i
	# mkdir ~/medvedevGroup/ecoli/errorRates$i
	cd ~/medvedevGroup/ecoli/errorRates$i
	# bash ~/medvedevGroup/bubblepopping/src/scriptToAutomateTasks.sh $i
	# bash ~/medvedevGroup/bubblepopping/src/scriptForValidation.sh
	bash ~/medvedevGroup/bubblepopping/src/scriptForAlign.sh
done
