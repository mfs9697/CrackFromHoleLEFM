function tests = test_SIF_auxiliary_fields
%TEST_SIF_AUXILIARY_FIELDS  Consistency tests for analytical LEFM auxiliary fields.
    tests = functiontests(localfunctions);
end


function testPlaneStrainConstitutiveConsistency(testCase)
    mat=make_mat(210e3,0.3,1);
    pts=[ ...
        1.0e-3,-2.4; ...
        4.0e-3,-1.1; ...
        1.0e-2,-0.2; ...
        2.0e-2, 0.4; ...
        5.0e-2, 1.3; ...
        8.0e-2, 2.7];

    modes=[1,0;0,1;1,0.7];

    worst=0;
    for im=1:size(modes,1)
        for k=1:size(pts,1)
            r=pts(k,1); th=pts(k,2);
            x1=r*cos(th); x2=r*sin(th);
            A=SIF_LEFM_auxiliary_fields(x1,x2,modes(im,1),modes(im,2),mat, ...
                'DerivativeMode','analytic');
            worst=max(worst,A.constitutive_mismatch);
        end
    end
    verifyLessThan(testCase,worst,1e-10);
end


function testAnalyticAgainstFiniteDifference(testCase)
    mat=make_mat(4e3,0.3,1);
    cases=[ ...
        3e-3,-2.1,1,0; ...
        7e-3,-0.7,0,1; ...
        2e-2, 0.3,1,0.4; ...
        5e-2, 2.2,0.8,-0.6];

    for k=1:size(cases,1)
        r=cases(k,1); th=cases(k,2); KI=cases(k,3); KII=cases(k,4);
        x1=r*cos(th); x2=r*sin(th);

        A=SIF_LEFM_auxiliary_fields(x1,x2,KI,KII,mat, ...
            'DerivativeMode','analytic');
        F=SIF_LEFM_auxiliary_fields(x1,x2,KI,KII,mat, ...
            'DerivativeMode','finite_difference','FDRelStep',1e-6,'FDAbsStep',1e-12);

        rel=norm(A.GradU-F.GradU,'fro')/max(norm(A.GradU,'fro'),eps);
        verifyLessThan(testCase,rel,2e-5);
    end
end


function testModeParity(testCase)
    mat=make_mat(210e3,0.3,1);
    r=0.01;
    th=0.9;

    IP=SIF_LEFM_auxiliary_fields(r*cos(th), r*sin(th),1,0,mat);
    IM=SIF_LEFM_auxiliary_fields(r*cos(th),-r*sin(th),1,0,mat);

    % Mode I: sigma11,sigma22 even; sigma12 odd.
    verifyEqual(testCase,IP.sig(1),IM.sig(1),'RelTol',1e-12);
    verifyEqual(testCase,IP.sig(2),IM.sig(2),'RelTol',1e-12);
    verifyEqual(testCase,IP.sig(3),-IM.sig(3),'RelTol',1e-12);

    IIP=SIF_LEFM_auxiliary_fields(r*cos(th), r*sin(th),0,1,mat);
    IIM=SIF_LEFM_auxiliary_fields(r*cos(th),-r*sin(th),0,1,mat);

    % Mode II: sigma11,sigma22 odd; sigma12 even.
    verifyEqual(testCase,IIP.sig(1),-IIM.sig(1),'RelTol',1e-12);
    verifyEqual(testCase,IIP.sig(2),-IIM.sig(2),'RelTol',1e-12);
    verifyEqual(testCase,IIP.sig(3), IIM.sig(3),'RelTol',1e-12);
end


function mat=make_mat(E,nu,ps)
    if ps==1
        c=E/((1+nu)*(1-2*nu));
        D=c*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
    else
        D=E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
    end
    mat=struct('E',E,'nu',nu,'Dmat',D,'D',D,'ps',ps);
end
