mpirun -np 1 ./bin/xhpcg datagen/ibmpg1.data > result/ibmpg1_1.txt
mpirun -np 2 ./bin/xhpcg datagen/ibmpg1.data > result/ibmpg1_2.txt
mpirun -np 4 ./bin/xhpcg datagen/ibmpg1.data > result/ibmpg1_4.txt
mpirun -np 6 ./bin/xhpcg datagen/ibmpg1.data > result/ibmpg1_6.txt
mpirun -np 8 ./bin/xhpcg datagen/ibmpg1.data > result/ibmpg1_8.txt
mpirun -np 10 ./bin/xhpcg datagen/ibmpg1.data > result/ibmpg1_10.txt
mpirun -np 12 ./bin/xhpcg datagen/ibmpg1.data > result/ibmpg1_12.txt

