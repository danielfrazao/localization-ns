#include <stdio.h>
#include <math.h>
#include <iostream>
#define WANT_STREAM 
#include "newmat10/newmatio.h"
#include "newmat10/newmat.h"

using namespace NEWMAT;

typedef struct {

	int Id_Pacote;	// Id do pacote enviado pelo no movel
	int Id_Fixo;  	// Id do no fixo
	
	double x_Fixo;	// coordenada x do no fixo 
	double y_Fixo;  // coordenada y do no fixo

	double x_Fixo_Calc; //coordenada x do fixo com erro de rssi 
	double y_Fixo_Calc; //coordenada y do fixo com erro de rssi
		
	double x_Movel;	// coordenada x do no movel
	double y_Movel; //coordenada y do no movel
	
	double dist;   		// distancia do no fixo que recebeu o pacote para o no movel 
	double dist_Calc;   	// distancia do no fixo que recebeu o pacote para o no movel Com erro de rssi 
	 
} Dados;

//função que lê as cordenadas e as transfere para uma Matriz

void uploadMatriz( Matrix& A, Matrix& A2, Dados *D, int num_pontos ){
	int i, j, n;
	double *c = new double[2*num_pontos]; //aloca vetor de coeficientes
	double *c2 = new double[2*num_pontos]; //aloca vetor de coeficientes
	
	n = num_pontos - 1;
	for( i = 0, j = 0; i < n; i++, j++ ) {
		c[j] = 2*( D[i].x_Fixo - D[n].x_Fixo );
		c2[j] = 2*( D[i].x_Fixo_Calc - D[n].x_Fixo_Calc );
		j++;
		c[j] = 2*( D[i].y_Fixo - D[n].y_Fixo );
		c2[j] = 2*( D[i].y_Fixo_Calc - D[n].y_Fixo_Calc );  
	}
	A << c;
	A2 << c2;

	delete c;  //libera memória
	delete c2;  //libera memória
}

void upload_b( ColumnVector& b, ColumnVector& b2, Dados *D, int num_pontos ){
	int i, n;
	double *v = new double[num_pontos];
	double *v2 = new double[num_pontos];
	n =  num_pontos - 1;
	for( i = 0; i < n; i++ ){
		v[i] = (pow(D[i].x_Fixo,2) - pow(D[n].x_Fixo,2)) + (pow(D[i].y_Fixo,2) - pow(D[n].y_Fixo,2)) + (pow(D[n].dist,2) - pow(D[i].dist,2)) ;
		v2[i] = (pow(D[i].x_Fixo_Calc,2) - pow(D[n].x_Fixo_Calc,2)) + (pow(D[i].y_Fixo_Calc,2) - pow(D[n].y_Fixo_Calc,2)) + (pow(D[n].dist_Calc,2) - pow(D[i].dist_Calc,2)) ;
	}
	b << v;
	b2 << v2;
	delete v;
	delete v2;
}

inline double distancia(double x,double y,double x2,double y2){

double distancia;
	
	distancia = (x2-x)*(x2-x)+(y2-y)*(y2-y);
	distancia = pow(distancia,0.5);
	return distancia;
}


int circle_circle_intersection(double x0, double y0, double r0,
                               double x1, double y1, double r1,
                               double *xi, double *yi,
                               double *xi_prime, double *yi_prime)
{
  double a, dx, dy, d, h, rx, ry;
  double x2, y2;

  /* dx e dy sao as distancias vertical e horizontal entre os centros dos circulos  */
  dx = x1 - x0;
  dy = y1 - y0;

  /* Determine the straight-line distance between the centers. */
  //d = sqrt((dy*dy) + (dx*dx));
  d = hypot(dx,dy); // Suggested by Keith Briggs

  /* Check for solvability. */
  if (d > (r0 + r1))
  {
    /* no solution. circles do not intersect. */
    return 0;
  }
  if (d < fabs(r0 - r1))
  {
    /* no solution. one circle is contained in the other */
    return 0;
  }

  /* 'point 2' is the point where the line through the circle
   * intersection points crosses the line between the circle
   * centers.  
   */

  /* Determine the distance from point 0 to point 2. */
  a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

  /* Determine the coordinates of point 2. */
  x2 = x0 + (dx * a/d);
  y2 = y0 + (dy * a/d);

  /* Determine the distance from point 2 to either of the
   * intersection points.
   */
  h = sqrt((r0*r0) - (a*a));

  /* Now determine the offsets of the intersection points from
   * point 2.
   */
  rx = -dy * (h/d);
  ry = dx * (h/d);

  /* Determina a interseção absoluta dos pontos. */
  *xi = x2 + rx;
  *xi_prime = x2 - rx;
  *yi = y2 + ry;
  *yi_prime = y2 - ry;

  return 1;
}


int PontoMaisProximoCirculo(double px, double py, double *x, double *y, double ray){
	
double lenght;
double x_aux = 0;
double y_aux = 0;


//{x: (px - x) / (length / ray) + x, y: (py - y) / (length / ray) + y}
	
	lenght =  (2 * 3.14 * ray);
	x_aux = *x;
	y_aux = *y;
	
	x_aux += (px - x_aux)/(lenght/ray); 
	y_aux += (py - y_aux)/(lenght/ray);
	
	*x = x_aux;
	*y = y_aux;
	
	
	return 1;
}


int primeiraLocalizacao(FILE *fp) {
	
#define TAM 10
double a = 0;
double b = 0;
int i,k;
double determinante;

int QtdD0 = 0;
int QtdD1 = 0;

int Id_Pacote_Aux;
int Id_Fixo_Aux;
double x_Fixo_Aux,y_Fixo_Aux;
double x_Fixo_Calc_Aux,y_Fixo_Calc_Aux;
double x_Movel_Aux,y_Movel_Aux;
double dist_Aux;
double dist_Aux2;


Dados D[TAM];

i = 0;

	//Leitura dos primeiros dados

	fscanf(fp,"%d",&D[i].Id_Pacote);
	fscanf(fp,"%d",&D[i].Id_Fixo);
	
	fscanf(fp,"%lf",&D[i].x_Fixo);
	fscanf(fp,"%lf",&D[i].y_Fixo);
	
	fscanf(fp,"%lf",&D[i].x_Fixo_Calc);
	fscanf(fp,"%lf",&D[i].y_Fixo_Calc);

	
	fscanf(fp,"%lf",&D[i].x_Movel);
	fscanf(fp,"%lf",&D[i].y_Movel);
	
	D[i].dist = distancia(D[i].x_Fixo,D[i].y_Fixo,D[i].x_Movel,D[i].y_Movel);
	D[i].dist_Calc = distancia(D[i].x_Fixo_Calc,D[i].y_Fixo_Calc,D[i].x_Movel,D[i].y_Movel);
	
	while (!feof(fp)){	
	// Faz a segunda leitura em variaveis auxiliares para comparacao posterior

	fscanf(fp,"%d",&Id_Pacote_Aux);
	fscanf(fp,"%d",&Id_Fixo_Aux);
	
	fscanf(fp,"%lf",&x_Fixo_Aux);
	fscanf(fp,"%lf",&y_Fixo_Aux);
	
	fscanf(fp,"%lf",&x_Fixo_Calc_Aux);
	fscanf(fp,"%lf",&y_Fixo_Calc_Aux);
	
	
	fscanf(fp,"%lf",&x_Movel_Aux);
	fscanf(fp,"%lf",&y_Movel_Aux);
	
	dist_Aux = distancia(x_Fixo_Aux,y_Fixo_Aux,x_Movel_Aux,y_Movel_Aux);
	dist_Aux2 = distancia(x_Fixo_Calc_Aux,y_Fixo_Calc_Aux,x_Movel_Aux,y_Movel_Aux);

		if(D[i].Id_Pacote == Id_Pacote_Aux){ // Se o pacote for o mesmo detectado pelo nó anterior , ele vai guardando em um vetor 
			
			i++;
			
			D[i].Id_Pacote = Id_Pacote_Aux;
			D[i].Id_Fixo = Id_Fixo_Aux;
			D[i].x_Fixo = x_Fixo_Aux;
			D[i].y_Fixo = y_Fixo_Aux;
			D[i].x_Fixo_Calc = x_Fixo_Calc_Aux;
			D[i].y_Fixo_Calc = y_Fixo_Calc_Aux;
			D[i].x_Movel = x_Movel_Aux;
			D[i].y_Movel = y_Movel_Aux;
			D[i].dist = dist_Aux;
			D[i].dist_Calc = dist_Aux2;
	
			Id_Pacote_Aux = -1;

		}//Fim da condição (D[i].Id_Pacote == Id_Pacote_Aux)

		if((D[i].Id_Fixo != Id_Pacote_Aux) && (Id_Pacote_Aux != -1)){ // caso o pacote do nó atual for diferente ao id do pacote da variável auxiliar  
			
			if ( i == 0){	
		
				printf("\n1 \n%d %d %lf %lf %lf %lf %lf %lf\n",D[i].Id_Pacote,D[i].Id_Fixo,D[i].x_Fixo,D[i].y_Fixo,D[i].x_Fixo_Calc,D[i].y_Fixo_Calc,D[i].x_Movel,D[i].y_Movel);

				D[i].Id_Pacote = Id_Pacote_Aux;
				D[i].Id_Fixo = Id_Fixo_Aux;
				D[i].x_Fixo = x_Fixo_Aux;
				D[i].y_Fixo = y_Fixo_Aux;
				D[i].x_Fixo_Calc = x_Fixo_Calc_Aux;
				D[i].y_Fixo_Calc = y_Fixo_Calc_Aux;
				D[i].x_Movel = x_Movel_Aux;
				D[i].y_Movel = y_Movel_Aux;
				D[i].dist = dist_Aux;
				D[i].dist_Calc = dist_Aux2;

			}		
			
			if ( i == 1 ) { 
				
			
				printf("\n2 \n%d %d %lf %lf %lf %lf %lf %lf\n",D[0].Id_Pacote,D[0].Id_Fixo,D[0].x_Fixo,D[0].y_Fixo,D[0].x_Fixo_Calc,D[0].y_Fixo_Calc,D[0].x_Movel,D[0].y_Movel);
				printf("%d %d %lf %lf %lf %lf %lf %lf\n\n",D[1].Id_Pacote,D[1].Id_Fixo,D[1].x_Fixo,D[1].y_Fixo,D[1].x_Fixo_Calc,D[1].y_Fixo_Calc,D[1].x_Movel,D[1].y_Movel);
				

				i = 0;

				D[i].Id_Pacote = Id_Pacote_Aux;
				D[i].Id_Fixo = Id_Fixo_Aux;
				D[i].x_Fixo = x_Fixo_Aux;
				D[i].y_Fixo = y_Fixo_Aux;
				D[i].x_Fixo_Calc = x_Fixo_Calc_Aux;
				D[i].y_Fixo_Calc = y_Fixo_Calc_Aux;
				D[i].x_Movel = x_Movel_Aux;
				D[i].y_Movel = y_Movel_Aux;
				D[i].dist = dist_Aux;
				D[i].dist_Calc = dist_Aux2;

			}
			
		
			if ( i >= 2 ) {

				i = 2;

				//printf(" \nO valor de i = %d\n",i);					
	
				Matrix A(i, 2);
				Matrix A2(i, 2);
				
				uploadMatriz( A, A2 , D , i+1 );
			
				//cout << "Matriz A\n";
				//cout << A;
				//printf("\n");
				
				determinante = A.Determinant();
			
					if (determinante != 0){
						ColumnVector b(i), y(i); //constrói os vetores b e y
						ColumnVector b2(i), y2(i); //constrói os vetores b e y
						
						upload_b( b, b2 , D, i+1 ); //carrega os valores de b;
					//	cout << "Vetor b\n";
					//	cout << b;
						y = A.i() * b; //y recebe a solução do sistema
						y2 = A2.i() * b2; //y recebe a solução do sistema
						
						printf("\n3 %d\n " , D[i].Id_Pacote); 
						cout << y;
						cout << y2;
						printf("%lf %lf\n",D[i].x_Movel,D[i].y_Movel);
						//cout << "Vetor solução:\n "<< y << endl;
						//printf("  \n %d Determinate calculado %.3lf \n\n",D[i].Id_Pacote,determinante); 
						QtdD1++;
					}
					
					if (determinante == 0){
					//	printf("Impossivel Calcular"); 
						QtdD0++;
					}

				i = 0 ;

				D[i].Id_Pacote = Id_Pacote_Aux;
				D[i].Id_Fixo = Id_Fixo_Aux;
				D[i].x_Fixo = x_Fixo_Aux;
				D[i].y_Fixo = y_Fixo_Aux;
				D[i].x_Fixo_Calc = x_Fixo_Calc_Aux;
				D[i].y_Fixo_Calc = y_Fixo_Calc_Aux;
				D[i].x_Movel = x_Movel_Aux;
				D[i].y_Movel = y_Movel_Aux;
				D[i].dist = dist_Aux;
				D[i].dist_Calc = dist_Aux2;
			}		


		}//Fim (D[i].Id_Fixo != Id_Pacote_Aux) && (Id_Pacote_Aux != -1)
		
		
	}//Fim While (!feof)  

}//Fim da Função que calcula a posição

int segundaLocalizacao(FILE *fp){

double determinante;

int status;
int i = 0;
int k = 0;

double x1_inter,y1_inter = 0;
double x2_inter,y2_inter = 0;

double x1_inter_Calc,y1_inter_Calc = 0;
double x2_inter_Calc,y2_inter_Calc = 0;

double d1,d2 = 0;
double d1_Calc,d2_Calc = 0;

Dados D[3];
Dados V[1];
Dados Aux;

	while (!feof(fp)){

		fscanf(fp,"%d",&status);
		
		if (status == 1){
		
			fscanf(fp,"%d",&Aux.Id_Pacote);
			fscanf(fp,"%d",&Aux.Id_Fixo);
			fscanf(fp,"%lf",&Aux.x_Fixo);
			fscanf(fp,"%lf",&Aux.y_Fixo);
			fscanf(fp,"%lf",&Aux.x_Fixo_Calc);
			fscanf(fp,"%lf",&Aux.y_Fixo_Calc);
			fscanf(fp,"%lf",&Aux.x_Movel);
			fscanf(fp,"%lf",&Aux.y_Movel);
			
			printf("\n %d \n %d %d %lf %lf %lf %lf %lf %lf  \n",status,Aux.Id_Pacote,Aux.Id_Fixo,Aux.x_Fixo,Aux.y_Fixo,Aux.x_Fixo_Calc,Aux.y_Fixo_Calc,Aux.x_Movel,Aux.y_Movel);
			
		}
		
		if (status == 2){

			
			for(i = 0 ; i < 2; i++){  
			
				fscanf(fp,"%d",&D[i].Id_Pacote);
				fscanf(fp,"%d",&D[i].Id_Fixo);
				fscanf(fp,"%lf",&D[i].x_Fixo);
				fscanf(fp,"%lf",&D[i].y_Fixo);
				fscanf(fp,"%lf",&D[i].x_Fixo_Calc);
				fscanf(fp,"%lf",&D[i].y_Fixo_Calc);
				fscanf(fp,"%lf",&D[i].x_Movel);
				fscanf(fp,"%lf",&D[i].y_Movel);
				D[i].dist = distancia(D[i].x_Fixo,D[i].y_Fixo,D[i].x_Movel,D[i].y_Movel);
				D[i].dist_Calc = distancia(D[i].x_Fixo_Calc,D[i].y_Fixo_Calc,D[i].x_Movel,D[i].y_Movel);
				
			
			}	
			
			circle_circle_intersection(D[0].x_Fixo, D[0].y_Fixo, 15, D[1].x_Fixo, D[1].y_Fixo, 15 ,&x1_inter, &y1_inter, &x2_inter, &y2_inter);
			circle_circle_intersection(D[0].x_Fixo_Calc, D[0].y_Fixo_Calc, 15, D[1].x_Fixo_Calc, D[1].y_Fixo_Calc, 15 ,&x1_inter_Calc, &y1_inter_Calc, &x2_inter_Calc, &y2_inter_Calc);
			//printf("Pontos de intersecao (%lf,%lf) - (%lf,%lf)",x1_inter,y1_inter,x2_inter,y2_inter);
			
			d1 = distancia(V[k].x_Fixo,Aux.y_Fixo,x1_inter,y1_inter);
			d2 = distancia(V[k].x_Fixo,Aux.y_Fixo,x2_inter,y2_inter);
			
			d1_Calc = distancia(V[k].x_Fixo_Calc,Aux.y_Fixo_Calc,x1_inter_Calc,y1_inter_Calc);
			d2_Calc = distancia(V[k].x_Fixo_Calc,Aux.y_Fixo_Calc,x2_inter_Calc,y2_inter_Calc);
			
			if(d1 < d2){
			
				D[i].x_Fixo = x1_inter; 
				D[i].y_Fixo = y1_inter;
				
			}
			
			if(d1_Calc < d2_Calc){
			
				D[i].x_Fixo_Calc = x1_inter_Calc; 
				D[i].y_Fixo_Calc = y1_inter_Calc;
				
			}
			
			if(d2 < d1){
			
				D[i].x_Fixo = x2_inter; 
				D[i].y_Fixo = y2_inter;
			}
			
			if(d2_Calc < d1_Calc){
			
				D[i].x_Fixo_Calc = x2_inter_Calc; 
				D[i].y_Fixo_Calc = y2_inter_Calc;
			}
			/** Garante que não vai haver pontos negativos **/
			if(D[i].x_Fixo_Calc < 0)
				D[i].x_Fixo_Calc = 0;
			if(D[i].y_Fixo_Calc < 0)
				D[i].y_Fixo_Calc = 0;
				
			printf("\n3\n%d %lf %lf %lf %lf %lf %lf\n ",D[0].Id_Pacote, D[i].x_Fixo, D[i].y_Fixo, D[i].x_Fixo_Calc, D[i].y_Fixo_Calc, D[0].x_Movel, D[0].y_Movel);
			
		}//Fim status == 2

		if (status == 3){
			
			fscanf(fp,"%d",&V[k].Id_Pacote);
			fscanf(fp,"%lf",&V[k].x_Fixo);
			fscanf(fp,"%lf",&V[k].y_Fixo);
			fscanf(fp,"%lf",&V[k].x_Fixo_Calc);
			fscanf(fp,"%lf",&V[k].y_Fixo_Calc);
			fscanf(fp,"%lf",&V[k].x_Movel);
			fscanf(fp,"%lf",&V[k].y_Movel);
			//Aux.dist = distancia(Aux.x_Fixo,Aux.y_Fixo,Aux.x_Movel,Aux.y_Movel);
			
			/** Garante que não vai haver pontos negativos **/
			if(V[k].x_Fixo_Calc < 0)
				V[k].x_Fixo_Calc = 0;
			if(V[k].y_Fixo_Calc < 0)
				V[k].y_Fixo_Calc = 0;
				
			printf("\n %d \n %d %lf %lf %lf %lf %lf %lf\n ",status, V[k].Id_Pacote, V[k].x_Fixo, V[k].y_Fixo, V[k].x_Fixo_Calc, V[k].y_Fixo_Calc, V[k].x_Movel, V[k].y_Movel);
			

		}


	}//Fim do While
	
//printf("\nDeterminante igual a zero : %d\n",QtdD0);
//printf("\nDeterminante diferente de zero : %d\n",QtdD1);
}

int TerceiraLocalizacao(FILE *fp){

double determinante;



int i = 0;

int status;

double x , y;
double x_Calc , y_Calc;
 
Dados D[3];
Dados V;
Dados Aux;

Aux.x_Fixo = 0;
Aux.y_Fixo = 0;

double dist_Aux;

	while (!feof(fp)){

		fscanf(fp,"%d",&status);
		
		if (status == 1){
		
			fscanf(fp,"%d",&Aux.Id_Pacote);
			fscanf(fp,"%d",&Aux.Id_Fixo);
			fscanf(fp,"%lf",&Aux.x_Fixo);
			fscanf(fp,"%lf",&Aux.y_Fixo);
			fscanf(fp,"%lf",&Aux.x_Fixo_Calc);
			fscanf(fp,"%lf",&Aux.y_Fixo_Calc);
			fscanf(fp,"%lf",&Aux.x_Movel);
			fscanf(fp,"%lf",&Aux.y_Movel);
			
			x = Aux.x_Fixo;
			y = Aux.y_Fixo;
			
			x_Calc = Aux.x_Fixo_Calc;
			y_Calc = Aux.y_Fixo_Calc;
			
			PontoMaisProximoCirculo(V.x_Fixo, V.y_Fixo, &x, &y, 15);
			PontoMaisProximoCirculo(V.x_Fixo_Calc, V.y_Fixo_Calc, &x_Calc, &y_Calc, 15);
			
			if(x_Calc < 0)
				x_Calc = 0;
			if(y_Calc < 0)
				y_Calc = 0;
				
			printf("\n 3 \n %d %lf %lf %lf %lf %lf %lf\n ", Aux.Id_Pacote, x, y,x_Calc, y_Calc, Aux.x_Movel, Aux.y_Movel);
			
		}
		
	
		
		if (status == 3){
			
			fscanf(fp,"%d",&V.Id_Pacote);
			fscanf(fp,"%lf",&V.x_Fixo);
			fscanf(fp,"%lf",&V.y_Fixo);
			fscanf(fp,"%lf",&V.x_Fixo_Calc);
			fscanf(fp,"%lf",&V.y_Fixo_Calc);
			fscanf(fp,"%lf",&V.x_Movel);
			fscanf(fp,"%lf",&V.y_Movel);
			//Aux.dist = distancia(Aux.x_Fixo,Aux.y_Fixo,Aux.x_Movel,Aux.y_Movel);
			
			if(V.x_Fixo_Calc < 0)
				V.x_Fixo_Calc = 0;
			if(V.y_Fixo_Calc < 0)
				V.y_Fixo_Calc = 0;
				
			printf("\n%d \n %d %lf %lf %lf %lf %lf %lf\n ",status, V.Id_Pacote, V.x_Fixo, V.y_Fixo, V.x_Fixo_Calc, V.y_Fixo_Calc, V.x_Movel, V.y_Movel);
			

		}

	}//Fim do While
	

}
int ResultadoLocalizacaoSegunda(FILE *fp){
int statusAux;
int Id_Pacote_Aux,Id_Fixo_Aux;
double x_Fixo_Aux,y_Fixo_Aux;
double x_Movel_Aux,y_Movel_Aux;
double x_Movel_Calc_Aux,y_Movel_Calc_Aux;
double x_Movel_Real, y_Movel_Real;
double MediaErro=0;
double MediaErro_Calc=0;
double dist_Aux;
double dist_Aux_Calc;
int num_loc = 0;
int descartados = 0;

	while (!feof(fp)){

		fscanf(fp,"%d",&statusAux);
			
			if (statusAux == 1){
				fscanf(fp,"%d",&Id_Pacote_Aux);
				fscanf(fp,"%d",&Id_Fixo_Aux);
				
				fscanf(fp,"%lf",&x_Movel_Aux);
				fscanf(fp,"%lf",&y_Movel_Aux);
				
				fscanf(fp,"%lf",&x_Movel_Calc_Aux);
				fscanf(fp,"%lf",&y_Movel_Calc_Aux);
				
				fscanf(fp,"%lf",&x_Movel_Real);
				fscanf(fp,"%lf",&y_Movel_Real);	
			}
		
			if(statusAux == 3){
			
				fscanf(fp,"%d",&Id_Pacote_Aux);
				fscanf(fp,"%lf",&x_Movel_Aux);
				fscanf(fp,"%lf",&y_Movel_Aux);
				
				fscanf(fp,"%lf",&x_Movel_Calc_Aux);
				fscanf(fp,"%lf",&y_Movel_Calc_Aux);
				
				fscanf(fp,"%lf",&x_Movel_Real);
				fscanf(fp,"%lf",&y_Movel_Real);			
				
				dist_Aux = distancia(x_Movel_Aux,y_Movel_Aux,x_Movel_Real,y_Movel_Real);
				dist_Aux_Calc = distancia(x_Movel_Calc_Aux,y_Movel_Calc_Aux,x_Movel_Real,y_Movel_Real);
				
				printf("\n%d %lf %lf \n",Id_Pacote_Aux,dist_Aux,dist_Aux_Calc);
				//if(dist_Aux<100){
				MediaErro = dist_Aux + MediaErro;
				MediaErro_Calc = dist_Aux_Calc + MediaErro_Calc;
				num_loc++;
				//}
				//else{
				//	descartados++;
			//	}
			}
	}
	MediaErro = (MediaErro/num_loc);
	MediaErro_Calc = (MediaErro_Calc/num_loc);
	
	printf("\nMedia de Erro : %lf \nMedia de erro com Rssi : %lf ",MediaErro,MediaErro_Calc);
	printf("\nNumero de pacotes detectados : %d",num_loc);
	
}

int ResultadoLocalizacao(FILE *fp){


int i = 0;
int IdMaior,IdMenor;
int IdMaior_Calc,IdMenor_Calc;
int status;
int Id_Pacote_Aux;
int Id_Fixo_Aux;
double x_Movel_Aux,y_Movel_Aux;
double x_Movel_Calc_Aux,y_Movel_Calc_Aux;
double x_Movel_Real, y_Movel_Real;

double MediaErro,MenorErro,MaiorErro,MediaErro_Calc,MaiorErro_Calc,MenorErro_Calc;
Dados D[3];
double dist_Aux;
double dist_Aux_Calc;

int checaID = -1;

MenorErro = 1000;
MaiorErro = -1000;

MenorErro_Calc = 1000;
MaiorErro_Calc = -1000;

	while (!feof(fp)){

		fscanf(fp,"%d",&status);

		if (status == 3){

			fscanf(fp,"%d",&Id_Pacote_Aux);
			fscanf(fp,"%lf",&x_Movel_Aux);
			fscanf(fp,"%lf",&y_Movel_Aux);
			fscanf(fp,"%lf",&x_Movel_Calc_Aux);
			fscanf(fp,"%lf",&y_Movel_Calc_Aux);
			fscanf(fp,"%lf",&x_Movel_Real);
			fscanf(fp,"%lf",&y_Movel_Real);			
			dist_Aux = distancia(x_Movel_Aux,y_Movel_Aux,x_Movel_Real,y_Movel_Real);
			dist_Aux_Calc = distancia(x_Movel_Calc_Aux,y_Movel_Calc_Aux,x_Movel_Real,y_Movel_Real);
			
			if(Id_Pacote_Aux != checaID ){
			
				printf("Id Pacote :%d \nx0: %lf y0: %lf x1: %lf y1: %lf \nx2: %lf y2: %lf \nErro : %lf\n\n",Id_Pacote_Aux, x_Movel_Aux,y_Movel_Aux,x_Movel_Calc_Aux,y_Movel_Calc_Aux,x_Movel_Real, y_Movel_Real, dist_Aux);
			
				if(dist_Aux > MaiorErro){
			
					MaiorErro = dist_Aux;
					IdMaior = Id_Pacote_Aux;
				}
				
				if(dist_Aux_Calc > MaiorErro_Calc){
				
					MaiorErro_Calc = dist_Aux_Calc;
					IdMaior_Calc = Id_Pacote_Aux;
				}  
				
				if(dist_Aux < MenorErro){
				
					MenorErro = dist_Aux;
					IdMenor = Id_Pacote_Aux;
				} 
			
				if(dist_Aux_Calc < MenorErro_Calc){
			
					MenorErro_Calc = dist_Aux_Calc;
					IdMenor_Calc = Id_Pacote_Aux;
				}
	

				MediaErro = dist_Aux + MediaErro;
				MediaErro_Calc = dist_Aux_Calc + MediaErro_Calc;
				i++;
				
				checaID = Id_Pacote_Aux;
			}
		}

	}//Fim do While
	
MediaErro = (MediaErro/i);

MediaErro_Calc = (MediaErro_Calc/i);

printf("\nNumero de Pacotes detectados : %d \n",i);

printf("\nErro Medio :%lf \nErro Medio com Rssi: %lf\n",MediaErro,MediaErro_Calc);

printf("\nMaior Erro :\nId : %d Erro: %lf \nErro com Rssi: %lf\n ",IdMaior,MaiorErro,MaiorErro_Calc);
printf("\nMenor Erro :\nId : %d Erro: %lf \nErro com Rssi: %lf\n ",IdMenor,MenorErro,MenorErro_Calc);

}

int main(int argc, char *argv[]){
	
	FILE *fp;   

	fp = fopen(argv[1],"r");
	freopen("PrimeiraLocalizacao","w",stdout);
	primeiraLocalizacao(fp);
	
	fp = fopen("PrimeiraLocalizacao","r");
	freopen("SegundaLocalizacao","w",stdout);
	segundaLocalizacao(fp);

	fp = fopen("SegundaLocalizacao","r");
	freopen("ResultadoSegunda","w",stdout);
	ResultadoLocalizacaoSegunda(fp);
	
	fp = fopen("SegundaLocalizacao","r");
	freopen("TerceiraLocalizacao","w",stdout);
	TerceiraLocalizacao(fp);

	fp = fopen("TerceiraLocalizacao","r");
	freopen("ResultadoLocalizacao","w",stdout);
	ResultadoLocalizacao(fp);

	return 0;		
	
}
