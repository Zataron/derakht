function advect()

    addpath src
    addpath src/semilag
    addpath src/tree
    addpath src/common

    % GLOBAL VARIABLES
    global VERBOSE;
    global GVFREQ;
    global DIM;
    global OM;

    global MAX_ERROR_PER_NODE;
    global MAX_LEVEL;
    global RES_PER_NODE;

    global INTERP_TYPE;
    global CHEB_IMPL;
    global VIS;

    global TINIT;
    global TN;
    global DT;

    global OUTPUT_DIR;

    VPREVTSTEP  = 1;
    VCURTSTEP   = 2;
    VNEXTSTEP   = 3;
    V2NEXTSTEP  = 4;

    % MAIN SCRIPT
    %fconc_exact   = @get_gaussian;
    fconc_exact   = @get_zalesak;
    %fconc_exact   = @get_cylinder_box;
    %fconc_exact   = @get_inviscid_burgers;
    fvel_exact    = @get_vel_exact;
    fvelx_exact   = @get_velx_exact;
    fvely_exact   = @get_vely_exact;
    %fdo_refine    = @do_refine;         % refinement criterion
    fdo_refine    = @do_refine_modified;         % refinement criterion

    % OUTPUT FOLDER
    timestamp = datestr(now,30);
    dir_name = [OUTPUT_DIR,'/res_',timestamp,'/'];
    mkdir(dir_name);

    diary([dir_name 'diary.txt']);

    t = TINIT + [-DT 0 DT 2*DT];

    % CONSTRUCT AND INIT THE TREES (VELOCITY AND CONCENTRAITION FEILDS)
    
    % VELOCITY (TIME-INDEPENDENT)
    u = qtree;
    u.insert_function(fvelx_exact,@do_refine,t(2));
    qdata.init_data(u,fvelx_exact,RES_PER_NODE,t(2));

    v = qtree;
    v.insert_function(fvely_exact,@do_refine,t(2));
    qdata.init_data(v,fvely_exact,RES_PER_NODE,t(2));

    % CONCENTRATION TREE
    cinit = qtree;
    cinit.insert_function(fconc_exact,@do_refine,TINIT);
    qdata.init_data(cinit,fconc_exact,RES_PER_NODE,TINIT);

    plot_tree(cinit,0); 
    
    % COMPUTE INITIAL TREES' ERRORS
    err = zeros(TN+1,2);
    err(1,1)= 0;
    err(1,2)= qdata.compute_error(cinit, fconc_exact, TINIT, INTERP_TYPE);

    % fprintf('DEPTH: %3d    Q: %3d   TN:%3d   DT:%8.2e\n', MAX_LEVEL, RES_PER_NODE, TN, DT);
    fprintf('INPUT ERROR: %12.2e\n', err(1,2));

    c = cinit;
    main_time = tic;
    for tstep =1:TN
        fprintf('======================================\n');

        % ONE STEP SL ADVECTION
        % ADVECT
        fprintf('-> Compute SL\n');
        ts_time = tic;
        cnext = sl_tree(c,u,v,t,fdo_refine,fconc_exact,fvel_exact);
        toc(ts_time)                

        % PLOT THE RESULT
        plot_tree(cnext,tstep);

        % COMPUTE THE ERROR
        e = qdata.compute_error(cnext, fconc_exact, t(VNEXTSTEP), INTERP_TYPE);
        err(tstep+1,1) = tstep;
        err(tstep+1,2) = e;
        format longE
        disp(['TN: ', num2str(tstep), '   Error: ', num2str(e,'%12.2e')]);
        
        % SAVE ERRORS AFTER FULL REVOLUTION
        %if mod(tstep,100) == 0
        %    % calculate and output mass ratio
        %    I = qdata.get_mass_ratio(cnext, fconc_exact, t(VCURTSTEP), INTERP_TYPE);
        %    fprintf('mass ratio: %f \n',I);     
        %    I2 = qdata.get_mass_ratio_squared(cnext, fconc_exact, t(VCURTSTEP), INTERP_TYPE);
        %    fprintf('mass ratio squared: %f \n',I2);                

            % calculate and output e_diss and e_disp
        %    [e1,e2,e3,e4] = qdata.get_interpolation_errors(cnext, fconc_exact, t(VCURTSTEP));
        %    fprintf('e_diss: %e \n',e1);
        %    fprintf('e_disp: %e \n',e2);              
            
        %    %revnum = sprintf('%i',tstep/100);
        %    %[X,Y] = qdata.grid_points(cnext);
        %    %[interp_values] = qdata.grid_data(cnext);    
        %    %real_values = fconc_exact(t(VCURTSTEP),X,Y,0);            
        %    %save(['testresults/test_results_cqmsl_rev',revnum,'.mat'],'e','e1','e2','I','I2','X','Y','interp_values','real_values');           
        %end        

        % PREPARE FOR THE NEXT STEP
        c = cnext;
        t = t + DT;       
    end % for time step
    tot_elapsed_time = toc(main_time);
    
    % calculate and output mass ratio
    I = qdata.get_mass_ratio(cnext, fconc_exact, t(VCURTSTEP), INTERP_TYPE);
    fprintf('mass ratio: %f \n',I);     
    I = qdata.get_mass(cnext, INTERP_TYPE);
    fprintf('mass of resulting tree: %f \n',I);
    I = qdata.get_mass(cinit, INTERP_TYPE);
    fprintf('mass of initial tree: %f \n',I);    
    
    % calculate and output e_diss and e_disp
    [e1,e2,e3,e4] = qdata.get_interpolation_errors(cnext, fconc_exact, t(VCURTSTEP));
    fprintf('e_diss: %e \n',e1);
    fprintf('e_disp: %e \n',e2);  
    fprintf('e_sum: %e \n',e3);  
    fprintf('e_total: %e \n',e4);  
    
    %plot cross-section    
    qdata.draw_cross_section(cnext,0.5,'cross_Y');
    
    % SAVE THE RESULTS
    frep = fopen([OUTPUT_DIR,'report.dat'],'a');
    % fprintf(frep,'%12s %5s %5s %10s %10s %10s %5s %12s %12s\n', ...
    %         'TOL', 'Q', 'MaxD', 'INTRP', 'CHEBIMPL', 'DT', 'TN', ...
    %         'InRLINF', 'OutRLINF');
    fprintf(frep,'%12.2e %5d %5d %10s %10s %10.2e %5d %12.2e %12.2e\n', ...
            MAX_ERROR_PER_NODE, RES_PER_NODE, MAX_LEVEL, INTERP_TYPE, CHEB_IMPL, DT, TN, ...
            err(1,2), err(end,2));

    % SAVE THE ERROR
    ferror = fopen([dir_name,'error.out'],'w');
    % fprintf(ferror,'# DEPTH: %3d    Q: %3d   TN:%3d   DT:%8.2e TOT_TIME: %f\n', ...
    %         MAX_LEVEL, RES_PER_NODE, TN, DT, tot_elapsed_time);
    fprintf(ferror,'%6s %12s\n','tn','LINF');
    fprintf(ferror,'%6d %12.2e\n',err');
    fclose(ferror);

    diary off;


    % PLOT THE RESULTS
    function plot_tree(cnext, tstep)
        if ~VIS, return; end;

        f = figure('Name','SEMI-LAG ADVECTION','visible','off');
        %f = figure('Name','SEMI-LAG ADVECTION');
        cnext.plottree(0.5);
        qdata.plot_data(cnext);
        colorbar;

        title(['tstep: ',num2str(tstep)]);
        b=sprintf('%4.3d',tstep);
        s_fig_name = [dir_name,'fig_',b]; % ,
        saveas(f,s_fig_name,'png');
        close(f);
    end

    %/* ************************************************** */
    function [val, funcval] = do_refine(node,func,t)
    % hRect = qtree.highlight_node(node, 2);
    % disp('pause:');
    % pause
    % delete(hRect);
        [val, funcval] = tree_do_refine(node,func,MAX_ERROR_PER_NODE,MAX_LEVEL,RES_PER_NODE,t);
    end

    %/* ************************************************** */
    function [val] = do_refine_modified(node,func,t)
    % hRect = qtree.highlight_node(node, 2);
    % pause;
    % delete(hRect);
        [val] = tree_do_refine_modified(node,func,MAX_ERROR_PER_NODE,MAX_LEVEL,RES_PER_NODE,t);
    end
end

%/* ************************************************** */
function value = get_gaussian(t,x,y,z)
    global OM;
    xc = 0.5;
    yc = 0.5;
    xci = 0.6;
    yci = 0.5;

    [alpha,RHO] = cart2pol(xci-xc,yci-xc);
    alphat = alpha + t*OM;

    [xct,yct] = pol2cart(alphat,RHO);
    xct = xct + xc;
    yct = yct + yc;

    theta = 0;
    sigmax = 0.06;
    sigmay = 0.06;
    value = gaussian(x,y,xct,yct,theta,sigmax,sigmay);
end

%/* ************************************************** */
function value = get_gaussian_2(t,x,y,z)
    global OM;
    xc = 0.5;
    yc = 0.5;
    theta = 0;
    sigmax = 0.05;
    sigmay = 0.05;

    xc1 = 0.75;
    yc1 = 0.75;
    [alpha1,RHO1] = cart2pol(xc1-xc,yc1-xc);
    alphat1 = alpha1 + t*OM;
    [xct,yct] = pol2cart(alphat1,RHO1);
    xct = xct + xc;
    yct = yct + yc;
    value1 = gaussian(x,y,xct,yct,theta,sigmax,sigmay);

    xc2 = 0.25;
    yc2 = 0.25;
    [alpha2,RHO2] = cart2pol(xc2-xc,yc2-xc);
    alphat2 = alpha2 + t*OM;
    [xct,yct] = pol2cart(alphat2,RHO2);
    xct = xct + xc;
    yct = yct + yc;
    value2 = gaussian(x,y,xct,yct,theta,sigmax,sigmay);

    value = value1 + value2;

end

%/* ************************************************** */
function value = get_zalesak(t,x,y,z)
    value = zalesak(x, y);
end

%/* ************************************************** */
function value = get_cylinder_box(t,x,y,z)
    xi = 0.0;
    xf = 1.0;
    domain_length = xf - xi;
    xcenter = (xf + xi) / 2;
    ycenter = xcenter;
    value = zeros(size(x));
    %box
    box_length = 0.5*domain_length;
    box_height = 2;    
    bmask = x <= xcenter+box_length/2 & x >= xcenter-box_length/2 & y <= ycenter+box_length/2 & y >= ycenter-box_length/2;
    value(bmask) = box_height;
    %cylinder
    cylinder_height = box_height + 2;
    cylinder_radius = 0.1*box_length;
    xd = x - xcenter;
    yd = y - ycenter;
    dist = sqrt(xd.*xd + yd.*yd);
    cmask = dist <= cylinder_radius;
    value(cmask) = cylinder_height;    
end

%/* ************************************************** */
function value = get_inviscid_burgers(t,x,y,z)
    value = 0.25 + 0.5*sin(0.5*pi*(x+y));
end

%/* ************************************************** */
function value = get_simple_function(t,x,y,z)
    %some simple function to test get_mass() correctness
    value = x.^2 + y.^2;
end

%/* ************************************************** */
function [u,v,w] = get_vel_exact(t,x,y,z)
    xc = 0.5;
    yc = 0.5;
    zc = 0.5;
    [u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
end

%/* ************************************************** */
function u = get_velx_exact(t,x,y,z)
    [u,v,w] = get_vel_exact(t,x,y,z);
end

%/* ************************************************** */
function v = get_vely_exact(t,x,y,z)
    [u,v,w] = get_vel_exact(t,x,y,z);
end
