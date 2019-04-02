#include "bitmap.cpp"
#include "fft.cpp"
#define OMP 0 //set true to output middle products

void Proc(Bitmap* a, double(*ftr)(double, double), const char *sfx){
	//initialize filter
	int n = 2 << (int)ceil(log2(a->height > a->width ? a->height : a->width));
	cmplx **x, **X;
	double **H;
	x = new cmplx*[n];
	X = new cmplx*[n];
	H = new double*[n];
	BMP256 *t;
	if(OMP) t = new BMP256(n, n);
	for(int i = 0; i < n; i++){
		x[i] = new cmplx[n];
		X[i] = new cmplx[n];
		H[i] = new double[n];
		int u = i > n/2 ? n - i : i;
		for(int j = 0; j < n; j++){
			int v = j > n/2 ? n - j : j;
			H[i][j] = ftr(u, v);
			
			if(i < a->height && j < a->width) x[i][j].re = a->px[i][j];
			else x[i][j].re = 0;
			
			if(OMP) t->px[i][j] = H[i][j] * 255.0;
		}
	}
	if(OMP){
		t->setName(sfx, "bmp");
		t->writeFile(n/2, n/2);
	}
	//2D FFT
	FFT fft(n);
	fft.FFT2(x, X);
	double M;
	if(OMP){
		M = 255.0 / log(X[0][0].mod());
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				int m = M * log(X[i][j].mod());
				if(m < 0) t->px[i][j] = 0;
				else t->px[i][j] = m;
			}
		}
		t->setName(a->name, "bmp");
		t->addNameSuffix("_f");	
		t->writeFile(n/2, n/2);
	}
	//frequency domain filtering
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			X[i][j].re *= H[i][j];
			X[i][j].im *= H[i][j];
			if(OMP){
				int m = M * log(X[i][j].mod());
				if(m < 0) t->px[i][j] = 0;
				else t->px[i][j] = m;
			}
		}
	}
	if(OMP){
		t->addNameSuffix(sfx);
		t->writeFile(n/2, n/2);	
	}
	//IFFT & recover
	fft.FFT2(X, x);
	BMP256 g(a->height, a->width);
	for(int i = 0; i < a->height; i++){
		for(int j = 0; j < a->width; j++){
			double m = x[n-i-1][n-j-1].re/n/n;
			if(m < 0) g.px[i][j] = 0;
			else if(m > 255) g.px[i][j] = 255;
			else g.px[i][j] = round(m);
		}
	}
	g.setName(a->name, "bmp");
	g.addNameSuffix(sfx);
	g.writeFile();	
}

