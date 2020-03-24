function dM = dualMT(t,M)

    global B1RMS
    global T_2B
    global Delta
    global g
    global R_A
    global M_0A
    global R
    global M_0B1
    global M_0B2
    global T_2A
    global R_B
    global T_D
    
    w_1 = 2*pi*4257.6*B1RMS;
    theta = 0:pi/2e3:pi/2;
    g = trapz(theta,sin(theta)*sqrt(2/pi).*(T_2B./abs(3*cos(theta).^2-1)).*...
        exp(-2*(Delta*T_2B./abs(3*cos(theta).^2-1)).^2));
    R_rfB = w_1*w_1*pi*g;

    dM = zeros(4,1);
    
    dM(1) = R_A*(M_0A-M(1)) - R*M_0B1*M(1) - R*M_0B2*M(1) + M_0A*R*M(2)...
        + M_0A*R*M(4) - ((w_1/Delta)^2)*M(1)/T_2A;
    
    dM(2) = R_B*(M_0B1-M(2)) - R*M_0A*M(2) + R*M_0B1*M(1) - R_rfB*M(2);
    
    dM(3) = - Delta*Delta*15*T_2B*T_2B*R_rfB*M(3) - M(3)/T_D;
    
    dM(4) = R_B*(M_0B2-M(4)) - R*M_0A*M(4) + R*M_0B2*M(1) - R_rfB*M(4);

end