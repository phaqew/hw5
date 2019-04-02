#include "bitmap.cpp"
#include "fft.cpp"
#define OMP 1 //set true to output middle products

void bwlp(Bitmap* a, double D0, int order = 2){
	int n = 2 << (int)ceil(log2(a->height > a->width ? a->height : a->width));
	cmplx **x, **X;
	double **H;
	x = new cmplx*[n];
	X = new cmplx*[n];
	H = new double*[n];
	BMP256 bw(n, n);
	for(int i = 0; i < n; i++){
		x[i] = new cmplx[n];
		X[i] = new cmplx[n];
		H[i] = new double[n];
		int u = i > n/2 ? n - i : i;
		for(int j = 0; j < n; j++){
			int v = j > n/2 ? n - j : j;
			if(i < a->height && j < a->width) x[i][j].re = a->px[i][j];
			else x[i][j].re = 0;
			double d = (u*u + v*v)/D0/D0, dd = 1;
			for(int k = 0; k < order; k++) dd *= d;
			H[i][j] = 1.0 / (1.0 + dd);
			bw.px[i][j] = H[i][j] * 255.0;
		}
	}
	bw.setName("Butterworth.bmp");
	char sfx[6] = "_00_0"; sfx[4]+=order; 
	sfx[1]+=(int)D0/10; sfx[2]+=(int)D0%10;
	bw.addNameSuffix(sfx);
	bw.writeFile(n/2, n/2);
	
	FFT fft(n);
	fft.FFT2(x, X);
	double M = 255.0 / log(X[0][0].mod());
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			int t = M * log(X[i][j].mod());
			if(t < 0) bw.px[i][j] = 0;
			else bw.px[i][j] = t;
		}
	}
	bw.setName(a->name, "bmp");
	bw.addNameSuffix("_f");	
	bw.writeFile(n/2, n/2);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			X[i][j].re *= H[i][j];
			X[i][j].im *= H[i][j];
			int t = M * log(X[i][j].mod());
			if(t < 0) bw.px[i][j] = 0;
			else bw.px[i][j] = t;
		}
	}
	bw.addNameSuffix("_bwlp");
	bw.addNameSuffix(sfx);
	bw.writeFile(n/2, n/2);	
	
	fft.FFT2(X, x);
	BMP256 g(a->height, a->width);
	for(int i = 0; i < a->height; i++){
		for(int j = 0; j < a->width; j++){
			double t = x[n-i-1][n-j-1].re/n/n;
			if(t < 0) g.px[i][j] = 0;
			else if(t > 255) g.px[i][j] = 255;
			else g.px[i][j] = round(t);
		}
	}
	g.setName(a->name, "bmp");
	g.addNameSuffix("_bwlp");
	g.addNameSuffix(sfx);
	g.writeFile();	
}

void gslp(Bitmap* a, double D0){
	int n = 2 << (int)ceil(log2(a->height > a->width ? a->height : a->width));
	cmplx **x, **X;
	double **H;
	x = new cmplx*[n];
	X = new cmplx*[n];
	H = new double*[n];
	BMP256 gs(n, n);
	for(int i = 0; i < n; i++){
		x[i] = new cmplx[n];
		X[i] = new cmplx[n];
		H[i] = new double[n];
		int u = i > n/2 ? n - i : i;
		for(int j = 0; j < n; j++){
			int v = j > n/2 ? n - j : j;
			if(i < a->height && j < a->width) x[i][j].re = a->px[i][j];
			else x[i][j].re = 0;
			double d = (u*u + v*v)/D0/D0;
			H[i][j] = exp(-0.5 * d);
			gs.px[i][j] = H[i][j] * 255.0;
		}
	}
	gs.setName("Gaussian.bmp");
	char sfx[4] = "_00";
	sfx[1]+=(int)D0/10; sfx[2]+=(int)D0%10;
	gs.addNameSuffix(sfx);
	gs.writeFile(n/2, n/2);
	
	FFT fft(n);
	fft.FFT2(x, X);
	double M = 255.0 / log(X[0][0].mod());
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			int t = M * log(X[i][j].mod());
			if(t < 0) gs.px[i][j] = 0;
			else gs.px[i][j] = t;
		}
	}
	gs.setName(a->name, "bmp");
	gs.addNameSuffix("_f");	
	gs.writeFile(n/2, n/2);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			X[i][j].re *= H[i][j];
			X[i][j].im *= H[i][j];
			int t = M * log(X[i][j].mod());
			if(t < 0) gs.px[i][j] = 0;
			else gs.px[i][j] = t;
		}
	}
	gs.addNameSuffix("_gslp");
	gs.addNameSuffix(sfx);
	gs.writeFile(n/2, n/2);	
	
	fft.FFT2(X, x);
	BMP256 g(a->height, a->width);
	for(int i = 0; i < a->height; i++){
		for(int j = 0; j < a->width; j++){
			double t = x[n-i-1][n-j-1].re/n/n;
			if(t < 0) g.px[i][j] = 0;
			else if(t > 255) g.px[i][j] = 255;
			else g.px[i][j] = round(t);
		}
	}
	g.setName(a->name, "bmp");
	g.addNameSuffix("_gslp");
	g.addNameSuffix(sfx);
	g.writeFile();	
}

int main(){
	PGM5 t1("test1.pgm");
	bwlp(&t1, 20, 1);
	bwlp(&t1, 20, 2);
	bwlp(&t1, 20, 3);
	bwlp(&t1, 50, 1);
	bwlp(&t1, 50, 2);
	bwlp(&t1, 50, 3);
	gslp(&t1, 20);
	gslp(&t1, 50);
	/*	
	TIF t2("test2.tif");
	Proc(&t2);
	PGM5 t3("test3.pgm");
	Proc(&t3);
	TIF t4("test4.tif");
	Proc(&t4);*/
}

