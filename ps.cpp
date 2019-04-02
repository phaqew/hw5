#include "bitmap.cpp"
#include "fft.cpp"

void Proc(Bitmap* a){
	int n = 2 << (int)ceil(log2(a->height > a->width ? a->height : a->width));
	cmplx **x, **X;
	x = new cmplx*[n];
	X = new cmplx*[n];
	for(int i = 0; i < n; i++){
		x[i] = new cmplx[n];
		X[i] = new cmplx[n];
		for(int j = 0; j < n; j++){
			if(i < a->height && j < a->width) x[i][j].re = a->px[i][j];
			else x[i][j].re = 0;
		}
	}
	FFT fft(n);
	fft.FFT2(x, X);
	double *p = new double[n/2];
	for(int i = 0; i < n/2; i++) p[i] = 0;
	for(int i = 0; i < n; i++){
		int x = i > n/2 ? n - i : i;
		for(int j = 0; j < n; j++){
			int y = j > n/2 ? n - j : j;
			int t = sqrt(x * x + y * y);
			if(t >= n/2) t = n/2 - 1;
			p[t] += X[i][j].mod() * X[i][j].mod();
		}
	}
	double M = 0;
	for(int i = 0; i < n/2 - 1; i++) if(M < p[i]) M = p[i];
	for(int i = 0; i < n/2; i++) p[i] *= (100.0 / M);
	BMP256 b(100, n/2);
	for(int i = 0; i < n/2; i++){
		for(int j = 0; j < 100 && j < p[i]; j++){
			b.px[j][i] = 255;
		}
	}
	for(int i = 0; ; i++)
		if((b.name[i] = a->name[i]) == '.'){
			b.name[i++] = '_'; b.name[i++] = 'p';
			b.name[i++] = 's'; b.name[i++] = '.';
			b.name[i++] = 'b'; b.name[i++] = 'm';
			b.name[i++] = 'p'; b.name[i] = 0; break;
		} //attach "_ps.bmp"
	b.writeFile();
	for(int i = 1; i < n/2; i++) p[i] += p[i-1];
	for(int i = 0; i < n/2; i++) p[i] /= p[n/2-1];
	char k[20];
	for(int i = 0; i < 20; i++)
		if((k[i] = a->name[i]) == '.'){
			k[i++] = '_'; k[i++] = 'p';
			k[i++] = 's'; k[i++] = '.';
			k[i++] = 'c'; k[i++] = 's';
			k[i++] = 'v'; k[i] = 0; break;
		} //attach "_ps.csv"
	std::ofstream ou(k, std::ios::out);
	for(int i = 0; i < n/2; i++) ou << i << "," << p[i] <<"\n";
	ou.close();
}

int main(){
	PGM5 t1("test1.pgm");
	Proc(&t1);
	TIF t2("test2.tif");
	Proc(&t2);
	PGM5 t3("test3.pgm");
	Proc(&t3);
	TIF t4("test4.tif");
	Proc(&t4);
}
