
#include <iostream>	//standart C++ kütüphanesi
#include <omp.h>	//openMP kütüphanesi
#include <ctime>	//tarih ve saat işlemleri için

#define N 1000	//matrix size
#define NUM_THREADS 8	//paralel thread sayısı

using namespace std;	//string veya vektör kullanımı için gerekli

//paralel matris çarpımı yapılır
void parallelCode(float* a, float* b, float* c) {

	int i, j, k;

	cout << "Threads calisma sirasi	==>	";

	omp_set_num_threads(NUM_THREADS);	//paralel çalışacak thread sayısı ayarlanır

	clock_t start = clock();	//programın başlangıç zaman bilgisi alınır

	//a, b, c matris verileri tüm bloklar tarafından erişilebilirken i, j, k verileri bloklara özeldir 
	#pragma omp parallel shared(a,b,c) private(i,j,k)	//paralel çalışacak olan kod bloğu 
	{
		int ID = omp_get_thread_num();	//threads numaraları alınır
		cout<<ID;

	#pragma omp for	schedule(static)	//Iterasyonlar varsayılan olarak eşit bir şekilde iş parçacıkları arasında paylaştırılır.
		for (i = 0; i < N; i = i + 1) {
			for (j = 0; j < N; j = j + 1) {
				//c[i * N + j] = 0;
				for (k = 0; k < N; k = k + 1) {
					c[i * N + j] = (c[i * N + j] + (a[i * N + k] * b[k * N + j]));
				}
			}
		}
	}

	double calismaSuresi = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;	//programın bitiş zaman bilgisi alınır
	cout << endl << "sure ==>	" << endl << calismaSuresi << endl;
}

//seri matris çarpımı yapar
void serialCode(float* a, float* b, float* c) {

	clock_t start = clock();

	int i, j, k;

	for (i = 0; i < N; i++) {

		for (j = 0; j < N; j++) {

			//c[i * N + j] = 0.0;

			for (k = 0; k < N; k++) {

				c[i * N + j] += a[i * N + k] * b[k * N + j];
			}
		}
	}
	
	double calismaSuresi = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
	cout << endl << "sure ==>	" << calismaSuresi << endl;
}

//çarpım ve sonuç matrislerine, çarpma işlemi öncesi değer ataması yapar
void fillingMatrix(float* a, float* b, float* c, float d) {

	int i, j, k;

	for (i = 0; i < N; i++) {

		for (j = 0; j < N; j++) {

			a[i * N + j] = d;
			b[i * N + j] = d;
			c[i * N + j] = 0;
		}
	}
}

void blockData(float* a, float* b, float* c) {

	clock_t start = clock();	//programın başlangıç zaman bilgisi alınır

	int i, j, k, ii, jj, kk;

	cout << "Threads calisma sirasi	==>	";

	omp_set_num_threads(NUM_THREADS);	//paralel çalışacak thread sayısı ayarlanır

	//a, b, c matris verileri tüm bloklar tarafından erişilebilirken i, j, k verileri bloklara özeldir 
	#pragma omp parallel shared(a,b,c) private(i,j,k,ii,jj,kk)	//paralel çalışacak olan kod bloğu 
	{
		int ID = omp_get_thread_num();	//threads numaraları alınır
		cout << ID;

	#pragma omp for	schedule(static)	//Iterasyonlar varsayılan olarak eşit bir şekilde iş parçacıkları arasında paylaştırılır.
		for (ii = 0; ii < N / 10; ++ii)
			for (jj = 0; jj < N / 10; ++jj)
				for (kk = 0; kk < N / 10; ++kk)
					for (i = 10 * ii; i < 10 * (ii + 1); ++i)
						for (j = 10 * jj; j < 10 * (jj + 1); ++j)
							for (k = 10 * kk; k < 10 * (kk + 1); ++k)
								c[i * N + j] += a[i * N + k] * b[k * N + j];
	}

	double calismaSuresi = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;	//programın bitiş zaman bilgisi alınır
	cout << endl << "sure ==>	" << endl << calismaSuresi << endl;
	cout << "--"<<c[86252]<<"--";
}

void yazdir(float* c) {

	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			cout << c[i * N + j] << " ";
		}
		cout << endl;
	}
}

int main()
{
	//dinamik matris tanımlamaları yapılır
	float* a = new float[N * N];
	float* b = new float[N * N];
	
	float* c = new float[N * N];	//çarpım sonucunun atanacağı sonuç matrisi tanımlanır

	float veri = 1.0;	//a ve b matrisinin barındıracağı veri

	fillingMatrix(a, b, c, veri);	//matrislere ilk değer atamaları yapılır

	blockData(a, b, c);

	//parallelCode(a, b, c);	//a ve b matrisi arasında openMP kullanılarak paralel çarpma işlemi yapılır

	//serialCode(a, b, c);	//a ve b matrisi arasında seri çarpma işlemi yapılır

	//istenildiği takdirde c yani çarpım sonucu yazdırılır fakat N nin çok büyük olduğu durumlarda düzgün bir görsel sonuç elde edilmez
	//yazdir(c);

	return 0;
	
}

