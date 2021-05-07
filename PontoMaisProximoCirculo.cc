/*
closestCirclePoint( px: Integer, py: Integer, x: Integer, y: Integer, ray: Double ): Object
Retorna um objeto contendo duas propriedades (x e y), que especificam o ponto limítrofe de um círculo em relação a um ponto.

px  coordenada x do ponto
py  coordenada y do ponto
x   coordenada x do ponto de origem do círculo
y   coordenada y do ponto de origem do círculo
ray raio do círculo (comprimento dividido por 2)

closestCirclePoint = function(px, py, x, y, ray){
	var tg = (x += ray, y += ray, 0);
	return function(x, y, x0, y0){return Math.sqrt((x -= x0) * x + (y -= y0) * y);}(px, py, x, y) > ray ?
		{x: Math.cos(tg = Math.atan2(py - y, px - x)) * ray + x, y: Math.sin(tg) * ray + y}
		//{x: (px - x) / (length / ray) + x, y: (py - y) / (length / ray) + y}
		: {x: px, y: py};
};



*/

#include <stdio.h>
#include <math.h>

int PontoMaisProximoCirculo(double px, double py, double x, double, double ray){
	
double lenght;
	
	lenght =  (2 * 3,14 * ray);
	*x = ((px - x)/(length/ray) + x);
	*y = ((py - y)/(length/ray) + y);
	
	return 1;
}

int main(){

}
