function Velocidade_Descolagem(a, b, c)

g=9.81;
S=0.985;
rho=1.2922;
miu=0.01;
C_l=0.75;

Z=zeros(20,19);

for j=1:20
    for i=5:24
    
        k=i-4;
        C_d=j*0.05;
        m=i/2;

        A_2=(1/m)*(a+0.5*rho*S*(miu*C_l-C_d));
        A_1=(1/m)*b;
        A_0=(1/m)*c-miu*g;

        % TODAS AS CONSTANTES NAO NULAS
        if A_0*A_1*A_2 ~= 0

            A = sqrt(abs(A_1*A_1-4*A_0*A_2))/(2*abs(A_2));
            B = -A_2*A;
            D = -A_1/(2*A_2);
            F = A/B;

            % CASO 1 - SOLUÇÃO TANGENTE HIPERBÓLICA
            if 4*A_0*A_2 < A_1*A_1
                C = atanh(-D/A);
                E = -F*log(cosh(C));

                u = @(t) A*tanh(B*t+C)+D;
                X = @(t) D*t+F*log(cosh(B*t+C))+E-55;
                T = fzero(X,10);
            end
            
            % CASO 2 - SOLUÇÃO TANGENTE
            if 4*A_0*A_2 > A_1*A_1
                C = atan(-D/A);
                E = (1/A_2)*log(abs(cos(C)));

                u = @(t) A*tan(A*A_2*t+C)+D;
                X = @(t) D*t+E-(1/A_2)*log(abs(cos(-B*t+C)))-55;
                T = fzero(X,10);
            end

            % CASO 3 - SOLUÇÃO LOGARÍTMICA
            if 4*A_0*A_2 == A_1*A_1
                c = -1/D;
                E = (1/A_2)*log(abs(1/D));

                u = @(t) 1/(A_2*t+C)+D;
                X = @(t) D*t+E-(1/A_2)*ln(abs(A_2*t+(1/D)))-55;
                T = fzero(X,10);
            end
        end

        % UMA CONSTANTE NULA
        if A_2 == 0 && A_0*A_1 ~= 0
            u = @(t) (A_0/A_1)*(e^(A_1*t)-1);
            X = @(t) (A_0/A_1)*((e^(A_1*t))/A_1-t)-A_0/((A_1)^2)-55;
            T = fzero(X,10);
        end
        
        if A_1 == 0 && A_0*A_2 ~= 0
            if A_0*A_2 > 0
                u = @(t) sqrt(A_0/A_2)*tan(sqrt(A_0/A_2)*A_2*t);
                X = @(t) (-1/A_2)*log(abs(cos(sqrt(A_0/A_2)*A_2*t)))-55;
                T = fzero(X,10);
            end
            if A_0*A_2 < 0
                u = @(t) sqrt(abs(A_0/A_2))*tanh(-A_2*sqrt(abs(A_0/A_2))*t);
                X = @(t) (-1/A_2)*log(cosh(-A_2*sqrt(abs(A_0/A_2))*t))-55;
                T = fzero(X,10);
            end
        end
        
        if A_0 == 0 && A_1*A_2 ~= 0
            printf('Não tem solução coerente com o problema');
            u(T) = 0;
        end
        
        % DUAS CONSTANTES NULAS
        if A_0 == 0 && A_1 == 0 && A_2 ~=0
            printf('Não tem solução coerente com o problema');
            u(T) = 0;
        end
        
        if A_0 == 0 && A_2 == 0 && A_1 ~=0
            printf('Não tem solução coerente com o problema');
            u(T) = 0;
        end
        
        if A_1 == 0 && A_2 == 0 && A_0 ~=0
            u = @(t) A_0*t;
            X = @(t) (A_0/2)/t^2-55;
            T = fzero(X,10);
        end
        
        % TRÊS CONSTANTES NULAS
        if A_0 == 0 && A_1 == 0 && A_2 == 0
           printf('Não tem solução coerente com o problema');
           u(T)=0;
        end
        
        Z(j,k) = u(T);      %Velocidade à descolagem
    end
end

X_mesh = 0.055:0.005:0.15;
Y_mesh = 2.5:0.5:12;

surf(X_mesh,Y_mesh,Z);
xlabel('C_d')
ylabel('m [kg]')
zlabel('u(55) [m/s]')