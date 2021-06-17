function a = plotter()

    a = 1;
    eqr = linspace(0.75,0.95,20);
    NOx = [];
    CO = [];
    
    valCO = 0.3*0.001; % kg/kg
    valNOx = 17.59*0.001; % kg/kg
    
    fuel = 0.855; % kg/s
    air = 34.06198; % kg/s     
    
    AFstoich = 14.79; % from Sofia
    
    ICAO_CO = [];
    ICAO_NOx = [];
    
    for i = 1:length(eqr)        
        [~,~,~,COf, NOxf] = reactor1('neo',3283120,1.1*805,eqr(i));
        CO(i) = COf;
        NOx(i) = NOxf;
        
        Air = (AFstoich/eqr(i))*fuel;
        ICAO_CO(i) = valCO*fuel/(Air+fuel);
        ICAO_NOx(i) = (valNOx*fuel/(Air+fuel))*1e6;      
    end  
    
    clf; %  clear figure
    
    COdif = (CO - ICAO_CO)./((ICAO_CO+CO)./2);
    COmax = max(COdif);
    
    NOxdif = (NOx - ICAO_NOx)./((ICAO_NOx+NOx)./2);
    NOxmax = max(NOxdif);
    
    subplot(3,1,1);
    hold on;
    plot(eqr,ICAO_CO,'b','LineWidth',1.5)
    plot(eqr,CO,'r','LineWidth',1.5)
    %plot(eqr,COdif,'r','LineWidth',1.5)
    legend('ICAO','Cantera')
    hold off;
    xlabel('Equivalence ratio (-)');
    ylabel('CO')
    
    subplot(3,1,2);
    hold on;
    %plot(eqr,ICAO_CO,'b','LineWidth',1.5)
    %plot(eqr,NOxdif,'r','LineWidth',1.5)
    plot(eqr,ICAO_NOx,'b','LineWidth',1.5)
    plot(eqr,NOx,'r','LineWidth',1.5)
    legend('ICAO','Cantera')
    hold off;
    xlabel('Equivalence ratio (-)');
    ylabel('NOx')
    
    subplot(3,1,3);
    hold on;
    %plot(eqr,ICAO_NOx,'b','LineWidth',1.5)
    %plot(eqr,NOx,'r','LineWidth',1.5)
    plot(eqr,abs(COdif),'r','LineWidth',1.5)
    plot(eqr, abs(NOxdif),'g','LineWidth',1.5)
    plot(eqr,abs(COdif)+abs(NOxdif),'b','LineWidth',1.5)
    legend('CO relative difference','NOx relative difference','Total relative difference')
    hold off;
    xlabel('Equivalence ratio (-)');
    ylabel('CO+NOx')
    
end