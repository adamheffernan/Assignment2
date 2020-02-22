function Assignment2_Part2(part)
%% Function for Part 2 of Assignment 2 4700. Adam Heffernan 100977570
% Set Part to 1 for 2A, 2 for 2B, 3 for 2C-1, 4 for 2C-2, and 5 for 2D

%% Part 2A
if part == 1
    %Variables
    clear
    clc
    len=30;
    width=20;%Height
    sigma_map=ones(width,len);
    F=zeros(1,len*width);
    G=zeros(len*width,len*width);
    sigma_map(((3/4)*(width)):width,(len/3):(2*(len/3)))=1e-2;
    sigma_map(1:(width/4),(len/3):(2*(len/3)))=1e-2;
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
                
                Nyu=j+1+(i-1)*width;
                Nxu=j+(i)*width;
                Nxd=j+(i-1-1)*width;
                
                Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                
                G(Updater,Updater)=-(Ty_Over+Tx_Over+Tx_Under);
                G(Updater,Nyu)=Ty_Over;
                G(Updater,Nxu)=Tx_Over;
                G(Updater,Nxd)=Tx_Under;
                
                
            elseif(j==width)%top
                
                Nyd=j-1+(i-1)*width;
                Nxu=j+(i)*width;
                Nxd=j+(i-1-1)*width;
                
                Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                
                G(Updater,Updater)=-(Ty_Under+Tx_Over+Tx_Under);
                G(Updater,Nyd)=Ty_Under;
                G(Updater,Nxu)=Tx_Over;
                G(Updater,Nxd)=Tx_Under;
                
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
    
    [Ex,Ey] = gradient(-V_2D);

    figure(4)
    contourf(V_2D,30)
    title('Voltage Contour Map')
    xlabel('X Coordinate')
    ylabel('Y Coordinate')
    c=colorbar;
    c.Label.String = 'Volts';
    
    figure(5)
    quiver(Ex,Ey)
    title('Electric Field Vector Map')
    xlabel('X Coordinate')
    ylabel('Y Coordinate')
    xlim([1 len])
    ylim([1 width])
    c=colorbar;
    c.Label.String = 'E (V/m^2)';
    
    figure(6)
    surf(sigma_map,'EdgeColor','none')
    title('Conductivity Map')
    xlabel('X Coordinate')
    ylabel('Y Coordinate')
    view(0,90)
    xlim([1 len])
    ylim([1 width])
    c=colorbar;
    c.Label.String = '\sigma (S/m)';
    
    figure(7)
    quiver(Ex.*sigma_map,Ey.*sigma_map)
    title('Current Density Map')
    xlabel('X Coordinate')
    ylabel('Y Coordinate')
    xlim([1 len])
    ylim([1 width])
    c=colorbar;
    c.Label.String = 'J_d (A/m^2)';

    %% Part 2B Multiple different mesh sizes 
    
elseif part == 2
        clear F
        clear G
        leftContact = (1:5);
        rightContact = (1:5);
        
    for multiplier=1:5 
       
        %Variables
        len=30*multiplier;
        width=20*multiplier;
        F=zeros(1,len*width);
        G=zeros(len*width,len*width);
        sigma_map=ones(width,len);
        sigma_map(((3/4)*(width)):width,(len/3):(2*(len/3)))=1e-2;
        sigma_map(1:(width/4),(len/3):(2*(len/3)))=1e-2;
        
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                
                if(i==1)%left
                    
                    G(updater,updater)=1;
                    F(updater)=1;
                    
                elseif(i==len)%right
                    
                    %G(N,:)=0;
                    G(updater,updater)=1;
                    F(updater)=0;
                    
                elseif(j==1)%bottom
                    Nyu=j+1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Over+Tx_Over+Tx_Under);
                    G(updater,Nyu)=Ty_Over;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                       
                elseif(j==width)%top
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Under+Tx_Over+Tx_Under);
                    G(updater,Nyd)=Ty_Under;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                    
                else%middle
                    Nyu=j+1+(i-1)*width;
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Over+Ty_Under+Tx_Over+Tx_Under);
                    G(updater,Nyu)=Ty_Over;
                    G(updater,Nyd)=Ty_Under;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                    
                end
            end
        end
        
        F = F';
        V=G\F;
        V_2D = zeros(width, len);
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                V_2D(j,i)=V(updater);
            end
        end
        
        [Ex,Ey]=gradient(-V_2D);
        Emag = sqrt(Ex.^2 + Ey.^2);
        leftContact(multiplier)=sum(sigma_map(:,1).*Emag(:,1));
        rightContact(multiplier)=sum(sigma_map(:,len).*Emag(:,len));
        
    end
    
    figure(8);
    hold on;
    plot(leftContact,'LineStyle','-','Marker','+','Color','black');
    plot(rightContact,'LineStyle','-.','Marker','o','Color','yellow');
    legend('Left Contact','Right Contact');
    hold on;
    title('Current v Mesh Size');
    ylabel('Relative Current');
    xlabel('Mesh size multiplier');

    
    %% Part 2C-1 Varrying Box Width
elseif part ==3
 %Variables 
    len=30;
    width=20;
    F=zeros(1,len*width);
    G=zeros(len*width,len*width);
    leftContact=zeros(1,50);
    rightContact=zeros(1,50);
    
    for varrying_width=1:width
        %Set conductivity map
        sigma_map=ones(width,len);
        sigma_map(1:varrying_width,(len/2):(2*(len/3)))=1e-2;
        sigma_map((width-varrying_width+1):width,(len/3):(2*(len/3)))=1e-2;
        
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                
                if(i==1)%left
                    
                    G(updater,updater)=1;
                    F(updater)=1;
                    
                elseif(i==len)%right
                    
                    G(updater,updater)=1;
                    F(updater)=0;
                    
                elseif(j==1)%bottom
                    Nyu=j+1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Over+Tx_Over+Tx_Under);
                    G(updater,Nyu)=Ty_Over;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                      
                elseif(j==width)%top
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Under+Tx_Over+Tx_Under);
                    G(updater,Nyd)=Ty_Under;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                    
                else%middle
                    Nyu=j+1+(i-1)*width;
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Over+Ty_Under+Tx_Over+Tx_Under);
                    G(updater,Nyu)=Ty_Over;
                    G(updater,Nyd)=Ty_Under;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                    
                end
            end
        end
        
        V=G\F';
        V_2D = zeros(width,len);
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;            
                V_2D(j,i)=V(updater);        
            end
        end
        
        [Ex,Ey]=gradient(-V_2D);   
        Emag = sqrt(Ex.^2 + Ey.^2);
        leftContact(varrying_width)=sum(sigma_map(:,1).*Emag(:,1));
        rightContact(varrying_width)=sum(sigma_map(:,len).*Emag(:,len));
        
    end
    figure(9)
    hold on;
    plot(leftContact,'LineStyle','-','Color','black');
    plot(rightContact,'LineStyle','-.','Color',[0.6350, 0.0780, 0.1840]);
    legend('Left Contact','Right Contact');
    xlim([0 20])
    hold off;
    title('Current v Varrying Width')
    xlabel('Width Variation from Original Mesh Size')
    ylabel('Relative Current')
    %% Part 2C-2 Varrying Box Height
elseif part == 4
    %Variables 
    len=30;
    width=20;
    F=zeros(1,len*width);
    G=zeros(len*width,len*width);
    leftContact =(1:len);
    rightContact =(1:len);
     for varrying_length=1:len
        %Set conductivity map
        sigma_map=ones(width,len);
        sigma_map(((3/4)*(width)):width,1:varrying_length)=1e-2;
        sigma_map(1:(width/4),(len-varrying_length+1):len)=1e-2;
        
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                
                if(i==1)%left
                    
                    G(updater,updater)=1;
                    F(updater)=1;
                    
                elseif(i==len)%right
                    
                    G(updater,updater)=1;
                    F(updater)=0;
                    
                elseif(j==1)%bottom
                    Nyu=j+1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Over+Tx_Over+Tx_Under);
                    G(updater,Nyu)=Ty_Over;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                      
                elseif(j==width)%top
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Under+Tx_Over+Tx_Under);
                    G(updater,Nyd)=Ty_Under;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                    
                else%middle
                    Nyu=j+1+(i-1)*width;
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Over+Ty_Under+Tx_Over+Tx_Under);
                    G(updater,Nyu)=Ty_Over;
                    G(updater,Nyd)=Ty_Under;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                    
                end
            end
        end
        
        V=G\F';
        V_2D = zeros(width,len);
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;            
                V_2D(j,i)=V(updater);        
            end
        end
        
        [Ex,Ey]=gradient(-V_2D);   
        Emag = sqrt(Ex.^2 + Ey.^2);      
        leftContact(varrying_length)=sum(sigma_map(:,1).*Emag(:,1));
        rightContact(varrying_length)=sum(sigma_map(:,len).*Emag(:,len));
        
    end
    figure(10)
    hold on;
    plot(leftContact,'LineStyle','-','Color','black');
    plot(rightContact,'LineStyle','--','Color',[0.9290, 0.6940, 0.1250]);
    legend('Left Contact','Right Contact');
    xlim([0 30])
    hold off;
    title('Current v Varrying Length')
    xlabel('Length Variation from Original Mesh Size')
    ylabel('Relative Current')
    

    %% Part 2D Varrying the conductivity density 
elseif part == 5
    %Variables 
    len=30;
    width=20;
    F=zeros(1,len*width);
    G=zeros(len*width,len*width);
    leftContact=zeros(1,50);
    rightContact=zeros(1,50);
    
    for multiplier=1:50
        %Set conductivity map
        sigma_map=ones(width,len);
        sigma_map(((3/4)*(width)):width,(len/3):(2*(len/3)))=1e-2*multiplier;
        sigma_map(1:(width/4),(len/3):(2*(len/3)))=1e-2*multiplier;
        
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                
                if(i==1)%left
                    
                    G(updater,updater)=1;
                    F(updater)=1;
                    
                elseif(i==len)%right
                    
                    G(updater,updater)=1;
                    F(updater)=0;
                    
                elseif(j==1)%bottom
                    Nyu=j+1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Over+Tx_Over+Tx_Under);
                    G(updater,Nyu)=Ty_Over;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                      
                elseif(j==width)%top
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Under+Tx_Over+Tx_Under);
                    G(updater,Nyd)=Ty_Under;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                    
                else%middle
                    Nyu=j+1+(i-1)*width;
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ty_Over=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Ty_Under=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Tx_Over=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Tx_Under=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ty_Over+Ty_Under+Tx_Over+Tx_Under);
                    G(updater,Nyu)=Ty_Over;
                    G(updater,Nyd)=Ty_Under;
                    G(updater,Nxu)=Tx_Over;
                    G(updater,Nxd)=Tx_Under;
                    
                end
            end
        end
        
        V=G\F';
        V_2D = zeros(width,len);
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;            
                V_2D(j,i)=V(updater);        
            end
        end
        
        [Ex,Ey]=gradient(-V_2D);   
        Emag = sqrt(Ex.^2 + Ey.^2);       
        leftContact(multiplier)=sum(sigma_map(:,1).*Emag(:,1));
        rightContact(multiplier)=sum(sigma_map(:,len).*Emag(:,len));
        
    end
    figure(11)
    hold on;
    plot(1e-2*(1:50),leftContact,'LineStyle','-','Color','black');
    plot(1e-2*(1:50),rightContact,'LineStyle','--','Color',[0, 0.75, 0.75]);
    legend('Left Contact','Right Contact');
    hold off;
    title('Current v Conductivity')
    xlabel('Conductivity')
    ylabel('Relative Current')
        
end