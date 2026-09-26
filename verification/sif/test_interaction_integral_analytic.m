function tests = test_interaction_integral_analytic
%TEST_INTERACTION_INTEGRAL_ANALYTIC
% Verify the interaction-integral algebra and normalization independently of FEM.
%
% We integrate two exact Williams K-fields over a polar annulus with a
% radial q-function.  For homogeneous isotropic elasticity the exact result
% is I = 2/Eeff * (KI*KauxI + KII*KauxII).

    tests=functiontests(localfunctions);
end


function testPureModeI(testCase)
    check_case(testCase,1.0,0.0,1.0,0.0);
end


function testPureModeII(testCase)
    check_case(testCase,0.0,1.0,0.0,1.0);
end


function testMixedModeSigns(testCase)
    check_case(testCase,1.30,-0.40,1.0,0.0);
    check_case(testCase,1.30,-0.40,0.0,1.0);
end


function check_case(testCase,KI,KII,KauxI,KauxII)
    E=4e3;
    nu=.3;
    mat=make_mat(E,nu);
    Eeff=E/(1-nu^2);

    rin=.01;
    rout=.03;
    nr=81;
    nth=801;

    r=linspace(rin,rout,nr);
    th=linspace(-pi+1e-6,pi-1e-6,nth);

    F=zeros(nr,nth);
    for ir=1:nr
        for it=1:nth
            x=r(ir)*cos(th(it));
            y=r(ir)*sin(th(it));

            A=SIF_LEFM_auxiliary_fields(x,y,KI,KII,mat, ...
                'DerivativeMode','analytic');
            B=SIF_LEFM_auxiliary_fields(x,y,KauxI,KauxII,mat, ...
                'DerivativeMode','analytic');

            qgrad=-(1/(rout-rin))*[cos(th(it));sin(th(it))];
            F(ir,it)=interaction_density(A,B,qgrad)*r(ir);
        end
    end

    I=trapz(r,trapz(th,F,2));
    Iexact=2/Eeff*(KI*KauxI+KII*KauxII);

    verifyEqual(testCase,I,Iexact,'RelTol',2e-5,'AbsTol',1e-10);

    Krecovered=.5*Eeff*I;
    Kexpected=KI*KauxI+KII*KauxII;
    verifyEqual(testCase,Krecovered,Kexpected,'RelTol',2e-5,'AbsTol',1e-8);
end


function d=interaction_density(actual,aux,qgrad)
    sig1=actual.sig;
    eps1=actual.eps;
    du1=actual.du_dx1;

    sig2=aux.sig;
    du2=aux.du_dx1;

    Wint=sig2.'*eps1;

    S1=[sig1(1),sig1(3);sig1(3),sig1(2)];
    S2=[sig2(1),sig2(3);sig2(3),sig2(2)];

    Avec=[ ...
        -Wint + dot(S1(:,1),du2) + dot(S2(:,1),du1); ...
                 dot(S1(:,2),du2) + dot(S2(:,2),du1)];

    d=Avec.'*qgrad;
end


function mat=make_mat(E,nu)
    c=E/((1+nu)*(1-2*nu));
    D=c*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
    mat=struct('E',E,'nu',nu,'Dmat',D,'D',D,'ps',1);
end
