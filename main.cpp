#include "filter.cpp"

double D0;
double idlpf(double u, double v){
	double d = u*u + v*v;
	if(d < D0*D0) return 1;
	else return 0;
}
double idhpf(double u, double v){
	double d = u*u + v*v;
	if(d > D0*D0) return 1;
	else return 0;
}
double midhpf(double u, double v){
	double d = u*u + v*v;
	if(d > D0*D0) return 2;
	else return 1;
}
int bworder;
double bwlpf(double u, double v){
	double d = (u*u + v*v)/D0/D0, dd = 1;
	for(int k = 0; k < bworder; k++) dd *= d;
	return 1.0 / (1.0 + dd);
}
double bwhpf(double u, double v){
	double d = (u*u + v*v)/D0/D0, dd = 1;
	for(int k = 0; k < bworder; k++) dd *= d;
	return dd / (1.0 + dd);
}
double mbwhpf(double u, double v){
	double d = (u*u + v*v)/D0/D0, dd = 1;
	for(int k = 0; k < bworder; k++) dd *= d;
	return dd / (1.0 + dd) + 1.0;
}
double gslpf(double u, double v){
	double d = (u*u + v*v)/D0/D0;
	return exp(-0.5 * d);
}
double gshpf(double u, double v){
	double d = (u*u + v*v)/D0/D0;
	return 1.0 - exp(-0.5 * d);
}
double mgshpf(double u, double v){
	double d = (u*u + v*v)/D0/D0;
	return 2.0 - exp(-0.5 * d);
}
double gain;
double lphpf(double u, double v){
	double d = u*u + v*v;
	return 4.0 * M_PI * M_PI * d / gain;
}
double mlphpf(double u, double v){
	double d = u*u + v*v;
	return 1.0 + 4.0 * M_PI * M_PI * d / gain;
}
double ldhpf(double u, double v){ //laplace with diagonal
	double d = u*u + v*v;
	return 12.0 * M_PI * M_PI * d / gain;
}
double mldhpf(double u, double v){ //laplace with diagonal
	double d = u*u + v*v;
	return 1.0 + 12.0 * M_PI * M_PI * d / gain;
}

int main(){
	PGM5 t1("test1.pgm");
	TIF t2("test2.tif");
	D0 = 20;
	Proc(&t1, idlpf, "_idlp_20");
	Proc(&t2, idlpf, "_idlp_20");
	bworder = 1;
	Proc(&t1, bwlpf, "_bwlp_20_1");
	Proc(&t2, bwlpf, "_bwlp_20_1");
	bworder = 2;
	Proc(&t1, bwlpf, "_bwlp_20_2");
	Proc(&t2, bwlpf, "_bwlp_20_2");
	bworder = 3;
	Proc(&t1, bwlpf, "_bwlp_20_3");
	Proc(&t2, bwlpf, "_bwlp_20_3");
	Proc(&t1, gslpf, "_gslp_20");
	Proc(&t2, gslpf, "_gslp_20");
	D0 = 50;
	Proc(&t1, idlpf, "_idlp_50");
	Proc(&t2, idlpf, "_idlp_50");
	bworder = 1;
	Proc(&t1, bwlpf, "_bwlp_50_1");
	Proc(&t2, bwlpf, "_bwlp_50_1");
	bworder = 2;
	Proc(&t1, bwlpf, "_bwlp_50_2");
	Proc(&t2, bwlpf, "_bwlp_50_2");
	bworder = 3;
	Proc(&t1, bwlpf, "_bwlp_50_3");
	Proc(&t2, bwlpf, "_bwlp_50_3");
	Proc(&t1, gslpf, "_gslp_50");
	Proc(&t2, gslpf, "_gslp_50");

	PGM5 t3("test3.pgm");
	TIF t4("test4.tif");
	D0 = 25;
	Proc(&t3, midhpf, "_midhp_25");
	Proc(&t4, midhpf, "_midhp_25");
	bworder = 1;
	Proc(&t3, mbwhpf, "_mbwhp_25_1");
	Proc(&t4, mbwhpf, "_mbwhp_25_1");
	bworder = 2;
	Proc(&t3, mbwhpf, "_mbwhp_25_2");
	Proc(&t4, mbwhpf, "_mbwhp_25_2");
	bworder = 3;
	Proc(&t3, mbwhpf, "_mbwhp_25_3");
	Proc(&t4, mbwhpf, "_mbwhp_25_3");
	Proc(&t3, mgshpf, "_mgshp_25");
	Proc(&t4, mgshpf, "_mgshp_25");
	D0 = 10;
	Proc(&t3, midhpf, "_midhp_10");
	Proc(&t4, midhpf, "_midhp_10");
	bworder = 1;
	Proc(&t3, mbwhpf, "_mbwhp_10_1");
	Proc(&t4, mbwhpf, "_mbwhp_10_1");
	bworder = 2;
	Proc(&t3, mbwhpf, "_mbwhp_10_2");
	Proc(&t4, mbwhpf, "_mbwhp_10_2");
	bworder = 3;
	Proc(&t3, mbwhpf, "_mbwhp_10_3");
	Proc(&t4, mbwhpf, "_mbwhp_10_3");
	Proc(&t3, mgshpf, "_mgshp_10");
	Proc(&t4, mgshpf, "_mgshp_10");
	gain = 512*512;
	Proc(&t3, mlphpf, "_mlphp");
	Proc(&t3, mldhpf, "_mldhp");
	gain = 1024*1024;
	Proc(&t4, mlphpf, "_mlphp");
	Proc(&t4, mldhpf, "_mldhp");
	return 0;
}

