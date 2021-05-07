all:
	g++ CalculaPosicao.cc -lm -lnewmat -o CalculaPosicao
	
	mv Saida-0.002 Saida
	./CalculaPosicao Saida
	mv ResultadoLocalizacao ResultadoLocalizacao-0.002
	mv Saida Saida-0.002
	
	mv Saida-0.003 Saida
	./CalculaPosicao Saida
	mv ResultadoLocalizacao ResultadoLocalizacao-0.003
	mv Saida Saida-0.003
	
	mv Saida-0.004 Saida
	./CalculaPosicao Saida
	mv ResultadoLocalizacao ResultadoLocalizacao-0.004
	mv Saida Saida-0.004
	
	mv Saida-0.005 Saida
	./CalculaPosicao Saida
	mv ResultadoLocalizacao ResultadoLocalizacao-0.005
	mv Saida Saida-0.005
	
	mv Saida-0.006 Saida
	./CalculaPosicao Saida
	mv ResultadoLocalizacao ResultadoLocalizacao-0.006
	mv Saida Saida-0.006
	
	mv Saida-0.007 Saida
	./CalculaPosicao Saida
	mv ResultadoLocalizacao ResultadoLocalizacao-0.007
	mv Saida Saida-0.007			
