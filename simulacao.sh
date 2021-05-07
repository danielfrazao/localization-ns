#!/bin/bash

echo "Digite o nome do diretório onde serão armazenados os resultados : "
read diretorio
echo "Digite o numero de simulações : "
read num_simulacoes

mkdir $diretorio

count=2;
aux=0;

while [ $aux != $num_simulacoes ]; do

	while [ $count != 10 ]; do
		
		./CalculaPosicao 0-Saida-0.00$count 
		
		mv ResultadoLocalizacao $aux-Resultado-0.00$count
		mv $aux-Resultado-0.00$count $diretorio/
		
		count=`expr $count + 1`;
	done
	
	count=2;
	
	aux=`expr $aux + 1`;
done



