function Assignment2_Part1(part)
%% Function for Part 1 of Assignment 2 4700. Adam Heffernan 100977570
%Input a 1 to simulate Part 1A or a 2 to run part 1B

%Variables
len=30;
width=20;
iter=100;
sigma_map=ones(width,len);
F=zeros(1,len*width);
G=zeros(len*width,len*width);
%% Part 1A
if part == 1
    for i=1:len
        for j=1:width
            Updater=j+(i-1)*width;
            
            if(i==1)%left
                
                G(Updater,Updater)=1;
                F(Updater)=1;
                
            elseif(i==len)%right
                
                G(Updater,Updater)=1;
                F(Updater)=0;
                
            elseif(j==1)%bottom
                
                G(Updater,Updater)=1;
                F(Updater)=0;
                
            elseif(j==width)%top
                
                G(Updater,Updater)=1;
                F(Updater)=0;
                
            else%middle
                Nyu=j+1+(i-1)*width;
                Nyd=j-1+(i-1)*width;
                Nxu=j+(i)*width;
                Nxd=j+(i-1-1)*width;
                
                Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                
                
                G(Updater,Updater)=-(Ty_Over+Ty_Under+Tx_Over+Tx_Under);
                G(Updater,Nyu)=Ty_Over;
                G(Updater,Nyd)=Ty_Under;
                G(Updater,Nxu)=Tx_Over;
                G(Updater,Nxd)=Tx_Under;
                
            end  
        end
    end
    
    F=F';
    
    V=G\F;
    V_1D = zeros(width,len);
    for i=1:len
        for j=1:width
            Updater=j+(i-1)*width;
            V_1D(j,i)=V(Updater);  
        end
    end
    
    figure(1)
    surf(V_1D,'EdgeColor','none','FaceColor','interp');
    title('Finite Difference Method for Single Boundary Conditon' );
    ylabel('Y');
    xlabel('X');
    c=colorbar;
    c.Label.String = 'Volts';
    view(0,90);
    
    %% Part 1 B
elseif part == 2
    for i=1:len
        for j=1:width
            Updater=j+(i-1)*width;
            
            if(i==1)%left
                
                G(Updater,Updater)=1;
                F(Updater)=1;
                
            elseif(i==len)%right
                
                G(Updater,Updater)=1;
                F(Updater)=1;
                
            elseif(j==1)%bottom
                
                G(Updater,Updater)=1;
                F(Updater)=0;
                
            elseif(j==width)%top
                
                G(Updater,Updater)=1;
                F(Updater)=0;
                
            else%middle
                Nyu=j+1+(i-1)*width;
                Nyd=j-1+(i-1)*width;
                Nxu=j+(i)*width;
                Nxd=j+(i-1-1)*width;
                
                Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
   
                G(Updater,Updater)=-(Ty_Over+Ty_Under+Tx_Over+Tx_Under);
                G(Updater,Nyu)=Ty_Over;
                G(Updater,Nyd)=Ty_Under;
                G(Updater,Nxu)=Tx_Over;
                G(Updater,Nxd)=Tx_Under;
                
            end
        end
    end
    
    F=F';
    
    V=G\F;
    V_2D = zeros(width,len);
    for i=1:len
        for j=1:width
            Updater=j+(i-1)*width;
            V_2D(j,i)=V(Updater);
        end
    end
    figure(2)
    surf(V_2D,'EdgeColor','none','FaceColor','interp');
    title('Finite Difference Method for Given Boundaries');
    ylabel('Y');
    xlabel('X');
    c=colorbar;
    c.Label.String = 'Volts';
    view(0,90);
    
    trunked_width=width-1;
    trunked_len=(len-1)/2.0;
    y_coordinates=0:1:trunked_width;
    x_coordinates=-trunked_len:1:trunked_len;
    [Xmesh,Ymesh]=meshgrid(x_coordinates,y_coordinates);
    V_analytical=zeros(length(y_coordinates),length(x_coordinates));

    for n=1:2:iter
        
        V_analytical=V_analytical+(4.*1./pi).*(cosh(n.*pi.*Xmesh./trunked_width).*sin(n.*pi.*Ymesh./trunked_width)./(n.*cosh(n.*pi.*trunked_len./trunked_width)));
        
    end
    figure(3)
    surf(V_analytical,'EdgeColor','none','FaceColor','interp')
    title('Analytical Solution of Laplaces Equation for Given Boundaries');
    ylabel('Y');
    xlabel('X');
    c=colorbar;
    c.Label.String = 'Volts';
    view(0,90);
else
    return
end
end