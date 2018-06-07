clc 
clear all 
close all
%vidObj=VideoWriter('flujo3.avi');           %Para definir y abrirr el video
%open(vidObj);

%% PAR�METROS DEL PROBLEMA  ***********************************************

nu=1.114e-6;      %nu del agua. 
l=1;              %Longitud de la cavidad en m. 
At=0.5;           %At inicial,En s
ro=9.99e2;        %Densidad del agua en Kg/m^3
vplaca=3.6/3.6;   %En m/s %A partir de unos 20 m/s se empieza a ir hacia la derecha con un paso de 0.5
iter=400;         %Iteraciones. N�mero de veces que se repite el bucle temporal independientemente del paso que demos
T0=10e-4;         %Tolerancia para el paso variable, determina el orden del error que queremos cometer en el esquema num�rico
                  %Esta toleracia debe ser siempre mayor que el ROUND OFF de la m�quina
Re=vplaca*l/nu
               
%% ************************************************************************
%% Para la extrapolacion de richardson esapcial  **************************

%for ni=1:3 
%ni=3
%if ni==1
%    N=5
%elseif ni==2
%    N=10
%elseif ni==3  
%    N=60
%end

%% ************************************************************************
%% DEFINI�N DE LA MALLA   *************************************************

N=40;            
n=N^2;             

for i=1:N  
    x(i)=(i-1)/N;    
end

y=x;
[X,Y]=meshgrid(x,y);

Ax=l/N;
Ay=Ax;
A=Ax;

%% ************************************************************************
%% DEFINICI�N DE LAS MATRICES DE CONDICIONES DE CONTORNO  *****************

ucont=zeros(N);          
ucont(1,:)=vplaca;       %Definici�n de la condici�n de velocidad superior

vcont=zeros(N);

ucont=sparse(ucont);     %Definici�n como sparse para ahorrar espacio (las sparse no guardan los ceros)
vcont=sparse(vcont);

%% ************************************************************************
%% DEFINICI�N DE LA MATRIZ DE VELOCIDADES (CONDICIONES INICIALES)

u = zeros(N);           %Condiciones de velocidad nula dentro dela cavidad
v = zeros(N);
u = sparse(u);          %Definici�n como sparse para ahorrar espacio
v = sparse(v);

%% ************************************************************************
%% DEFINICI�N DE LAS MATRICES DE DIFERENCICACI�N

%Definici�n de las Matrices D1 y D2 para la derivaci�n 
e = ones(N,1);
D1 = spdiags([-e 0*e e], -1:1, N, N);                  %Matriz de diferenciaci�n D1 de la forma: -1 0 1 
D2 = spdiags([e -2*e e], -1:1, N, N);                  %Matriz de diferenciaci�n D2 de la forma: 1,-2,1 

%%Definici�n de la matriz L para la ecuacci�n de poisson
e2 = ones(n,1);
L = (spdiags([e2 -4*e2 e2], -1:1, n, n))+...             %puntos suspensivos para que no se vaya la linea mucho hacia la derecha
(spdiags([e2], -N, n, n))+(spdiags([e2], N, n, n));      %Matriz de diferenciaci�n L para la ecuaci�n de poisson.
Dx = (spdiags([e2], -N, n, n))+(spdiags([e2], N, n, n)); %Matrix de diferenciaci�n Dx de N^2xN^2 (para la estabilidad)
Dy = spdiags([e2 0*e2 -e2], -1:1, n, n);                 %Matrix de diferenciaci�n Dx de N^2xN^2 (para la estabilidad)

for i=1:N-1 
    L(N*i,N*i+1)=0;                                    %Para terminar de definir L
end

%% ************************************************************************
%% BUCLE TEMPORAL

for f=1:iter
%% DEFINICI�N DE LA MATRIZ W DE VORTICIDAD EN EL PASO N A PARTIR DE LA VELOCIDAD U Y V EN EL PASO N

w=(1/(2*Ax))*(v*D1'-vcont)-(1/(2*Ay))*(D1*u-ucont);

%% ************************************************************************
%% CONDICIONES DE CONTORNO DE LA VORTICIDAD:  *****************************
%% VECTORES DE CONTORNO
    
wcontW(1)=0;
wcontW(N)=0;

wcontE(1)=0;
wcontE(N)=0;

wcontN(1)=0;
wcontN(N)=0;

wcontS(1)=0;
wcontS(N)=0;

for i=2:N-1
   wcontW(i)=(1/(2*Ax))*(-3*vcont(i,1)+4*v(i,1)-v(i,2))-(1/(2*Ay))*(ucont(i+1,1)-ucont(i-1,1));
   wcontE(i)=(1/(2*Ax))*(v(i,N-1)-4*v(i,N)+3*vcont(i,N))-(1/(2*Ay))*(ucont(i+1,N)-ucont(i-1,N));
end 
for j=2:N-1
   wcontN(j)=(1/(2*Ax))*(vcont(1,j+1)-vcont(1,j-1))-(1/(2*Ay))*(-3*ucont(1,j)+4*u(i,j)-u(1,j+1));
   wcontS(j)=(1/(2*Ax))*(vcont(N,j+1)-vcont(N,j-1))-(1/(2*Ay))*(u(N,j)-4*u(N,j)+3*ucont(N,j));
end
    
%%%%PARA FINALMENTE DEFINIR LAS MATRICES DE CONTORNOS
    
wcont1WE(N,N)=0;
wcont1WE(:,1)=-wcontW;
wcont1WE(:,N)=wcontE;

wcont1NS(N,N)=0;
wcont1NS(1,:)=-wcontN;
wcont1NS(N,:)=wcontS;

wcont2WE(N,N)=0;
wcont2WE(:,1)=wcontW;
wcont2WE(:,N)=wcontE;

wcont2NS(N,N)=0;
wcont2NS(1,:)=wcontN;
wcont2NS(N,:)=wcontS;

%% ***********************************************************************
%% PARA LA EXTRAPOLACI�N ESPACIAL DE RICHARSON ***************************
    
    %%EXTRAPOLACION TEMPORAL DE RICHARDSON
     %if t==100*At
     %    if N==5
     %        Fr1=-(1/2*Ax)*((w*D1'+wcont1WE).*u)-(1/2*Ay)*((D1*w+wcont1NS).*v)+...(1/Re)*(((1/Ax^2)*w*D2+wcont2WE)+((1/Ay^2)*D2*w+wcont2NS));
     %        FR1=Fr1(:);
     %    elseif N==10
     %        Fr2=-(1/2*Ax)*((w*D1'+wcont1WE).*u)-(1/2*Ay)*((D1*w+wcont1NS).*v)+(1/Re)*(((1/Ax^2)*w*D2+wcont2WE)+((1/Ay^2)*D2*w+wcont2NS));
     %        FR2=Fr2(:);
     %    elseif N==20
     %        Fr3=-(1/2*Ax)*((w*D1'+wcont1WE).*u)-(1/2*Ay)*((D1*w+wcont1NS).*v)+(1/Re)*(((1/Ax^2)*w*D2+wcont2WE)+((1/Ay^2)*D2*w+wcont2NS));
     %        FR3=Fr3(:);
     %    elseif N==40
     %        Fr4=-(1/2*Ax)*((w*D1'+wcont1WE).*u)-(1/2*Ay)*((D1*w+wcont1NS).*v)+(1/Re)*(((1/Ax^2)*w*D2+wcont2WE)+((1/Ay^2)*D2*w+wcont2NS));
     %        FR4=Fr4(:);
     %    end
     %end

%%  ***********************************************************************
%% C�LCULO DE LA F  *******************************************************
     
      F=-(1/2*Ax)*((w*D1'+wcont1WE).*u)-(1/2*Ay)*((D1*w+wcont1NS).*v)+(1/Re)*(((1/Ax^2)*w*D2+wcont2WE)+((1/Ay^2)*D2*w+wcont2NS));
      
      F=sparse(F);                    %sparse para seguir ahorrando espacio
      

%%  ***********************************************************************
%% ESQUEMA TEMPORAL AB3 PREDICTOR *****************************************
%%Para guardar espacio, a cada iteraci�n se va borrando Faux(1) con lo que Faux(2) pasa a ser 
%%Faux(1). (VER L�nea 206). Por eso se va llamando a la Faux(:,2) y Faux(:,3)y no a la Faux(:,f) y Faux(:,f-1)

    if (f==1) 
       Faux(:,1)=F(:);     
       w=w+At*F;                                      %Empezamos con Euler
    elseif (f==2)
       Faux(:,2)=F(:);
       w(:)=w(:)+(At/2)*[3*Faux(:,2)-Faux(:,1)];           %Adam basford 2
    else
       Faux(:,3)=F(:);
       w(:)=w(:)+(At/12)*[23*Faux(:,3)-16*Faux(:,2)+5*Faux(:,1)];     %AB3
    end

%% ************************************************************************
%% ESQUEMA TEMPORAL AM2 CORRECTOR *****************************************

    F=-(1/2*Ax)*((w*D1'+wcont1WE).*u)-(1/2*Ay)*((D1*w+wcont1NS).*v)+(1/Re)*(((1/Ax^2)*w*D2+wcont2WE)+((1/Ay^2)*D2*w+wcont2NS));
    if (f==1)
        wt=w;                            %En el primer paso no se puede hacer todav�a la correcci�n,da igual pq el error del arraque termina desapareciendo
    elseif (f==2)
        wt(:)=w(:)+(At/2)*(Faux(:,2)+F(:));                %AM1 
    elseif (f>2)
        wt(:)=w(:)+(At/12)*(8*Faux(:,3)+5*F(:)-Faux(:,2)); %AM2
        Faux(:,1)=[];                                       %En vez de guardar toda Faux, vamos quitando la primera columna para guardar espacio
    end
    
%% ESTIMACI�N DEL ERROR LOCAL:TEOREMA DE MILNE ****************************
      
      Ee=((-1/24)/((3/8)+(1/24)))*(w-wt);   %Estimaci�n del error 
      Ee=full(Ee);                          %de sparse a entera para poder hallar la norma
      w=wt+Ee;                               %Mejoramos la soluci�n a�adiendole laestimaci�n del error
      Eem=((-1/24)/((3/8)+(1/24)))*(w-wt);  %Estimaci�n del error entre la soluci�n mejorada y sin mejorar (Para comparar la mejora)
      E(f)=norm(Ee);                        %Norma de la estimaci�n del error local para poder representarlo
      E1(f)=norm(Eem);                      %Norma de la estimaci�n del error con la soluci�n mejorada.

%% PASO VARIABLE **********************************************************
%%Hemos definido una tolerancia para el error local T0. Para esa tolerancia
%%hay un At m�ximo asociado. Si queremos obtener esa tolerancia habr� que
%%calcular ese Atm�x y aplic�rselo a nuestro problema
%%Puesto en verde para no hacernos lios

      %if E<T0
      %  Atmax=(T0/E(f))^(1/3)*At;
      %   At=Atmax;
      %else
      %    Atmax=(T0/E(f))^(1/3)*At;
      %    At=Atmax;
      %end
      
%%  ***********************************************************************
%% SE CALCULA LA FUNCI�N DE CORRIENTE EN EL PASO N+1:
    
    psi=zeros(N);
    psi(:)=L\(-w(:)*A^2);
    
%% ************************************************************************
%% SE CALCULAN LAS VELOCIDADES U Y V EN EL PASO N+1:
    
    u=(1/(2*A))*D1*psi;
    v=(-1/(2*A))*psi*D1';

%% ************************************************************************
%% C�LCULO DEL M�DULO DE LA VELOCIDAD EN CADA PUNTO ***********************

    for i=1:N
        for j=1:N
        U(i,j)=sqrt(u(i,j)^2+v(i,j)^2);
       end
    end

%% ************************************************************************
%% REPRESENTACI�N DE LA SOLUCI�N  *****************************************

%if t==p                %poner solo si se usa la extrapolaci�n de richardon
%fig1 = figure(1); 
%delete(gca) 
%figure(1)
%streamslice(x,y,u,v);                               %Lineas de corriente
%set(fig1,'DoubleBuffer','on')
%set(gca,'nextplot','replace','Visible','on','Ydir','Reverse')
%figure(2)
%pcolor(x,y,U)                                      %M�dulo de la velocidad
%set(fig1,'DoubleBuffer','on')
%set(gca,'nextplot','replace','Visible','on','Ydir','Reverse')
%frame = getframe(gcf);                             
%writeVideo(vidObj,frame);                          %Para hacer el v�deo
end

%end
%close(vidObj);

%% ************************************************************************
%% FIGURA DE LA REPRESENTACI�N DE LA ESTIMACI�N ERROR *********************

%figure(4)
%hold on
%plot(E)
%plot(E1,'r')
%xlabel('Iteracione')
%ylabel('Estimaci�n del error')
%title('Estimaci�n del error local(M�todo de milne)')
%legend('Estimaci�n sin mejorar','Estimaci�n mejorada')

%%  ***********************************************************************
%% ERROR RICHARDSON ESPACIAL REPRESENTACI�N  ******************************

%for i=1:1:5
%     T(i)=(FR2(1+(2*(i-1)))-FR1(i))/(1-(1/4)) ;
%end
%Trichardson(1)=norm(T);
%for i=1:1:10
%      T(i)=(FR3(1+(2*(i-1)))-FR2(i))/(1-(1/4));
%end
%Trichardson(2)=norm(T);
%for i=1:1:20
%      T(i)=(FR4(1+(2*(i-1)))-FR3(i))/(1-(1/4));
%end
%Trichardson(3)=norm(T);
%N_vector=[5,10];

%figure(4)
%plot(N_vector,Trichardson)
%xlabel('N')
%ylabel('Trichadson')
%title('Error espacial')

%%  ***********************************************************************
%% CALCULO DE LA PRESI�N  *************************************************

lap=((u*D1')/(2*A))^2+((u.*(u*D2))/(A^2))+((v*D1')/(2*A))^2+((v.*(v*D2))/A^2)+...
((D1*u)/(2*A))^2+((u.*(D2*u))/(A^2))+((D1*v)/(2*A))^2+((v.*(D2*v))/A^2);

p=zeros(N);
p(:)=L\(-lap(:)*A^2);
[Fx,Fy]=gradient(p);                   %C�lculo del gradiente de la presi�n

%% REPRESENTACI�N DE LA PRESI�N  ******************************************
%%Mapa de intensidades de la presi�n. Recomendable acompa�arlo viendo el
%%gradiente de la presi�n

figure(5)
pcolor(x,y,p)
set(gca,'YDir','Reverse')
%set(h,'EdgeColor', 'none')

%% REPRESENTACI�N DEL GRADIENTE DE LA PRESI�N *****************************
%%con esta figura se puede ver la DIRECCI�N en la que act�a la presi�n

figure(6)
quiver(x,y,Fx,Fy)       
set(gca,'YDir','Reverse')
