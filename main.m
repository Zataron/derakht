function main()
    clear all;
    clear globals; % constants  and preamble

    % GLOBAL VARIABLES
    global VERBOSE;
    global GVFREQ;
    global DIM;
    global OM;
    global VIS;

    global MAX_ERROR_PER_NODE;
    global MAX_LEVEL;
    global RES_PER_NODE;
    global ADAPTIVE;

    global INTERP_TYPE;
    global CHEB_IMPL;
    global CHEB_KIND;
    global ENABLE_QMSL
    global ENABLE_CQMSL

    global TINIT;
    global TN;
    global DT;

    global OUTPUT_DIR;
    global TESTNUM;     %TODO: remove later on
    global FILTER;      %TODO: remove later on

    %ugrid_temporal_conv()
    ugrid_spatial_temporal_conv()

    function ugrid_temporal_conv()
    % OUTPUT FOLDER
        timestamp  = datestr(now,30);
        OUTPUT_DIR = ['./results/','ugrid_temporal_conv', '_', timestamp,'/'];
        mkdir(OUTPUT_DIR);

        frep = fopen([OUTPUT_DIR,'report.dat'],'a');
        fprintf(frep,'%12s %5s %5s %10s %10s %10s %5s %12s %12s\n', ...
                'TOL', 'Q', 'MaxD', 'INTRP', 'CHEBIMPL', 'DT', 'TN', ...
                'InRLINF', 'OutRLINF');

        % SIMULATION PARAMETERS
        VERBOSE = false;
        GVFREQ  = 0;
        OM      = 1;
        DIM     = 2;
        VIS     = true;

        % SPATIOAL RESOLUTION
        MAX_ERROR_PER_NODE = 1e-30; % Error per box
        MAX_LEVEL          = 4;
        RES_PER_NODE       = 14;    % Resolution per Node
        ADAPTIVE           = false;
        INTERP_TYPE        = 'cubic';
        %INTERP_TYPE        = 'CHEBYSHEV';
        %CHEB_IMPL          = 'CHEBFUN';
        %CHEB_IMPL          = 'IAS';
        %CHEB_KIND          = 2;

        % TEMPORAL RESOLUTION
        TINIT   = 0;
        tfinal  = pi/4;
        tn_list = 2.^[1:2];

        for lvl =1:length(tn_list)
            TN      = tn_list(lvl);
            DT      = tfinal/TN;
            % RUN
            advect();
        end
    end

    function ugrid_spatial_temporal_conv()
    % OUTPUT FOLDER
        timestamp  = datestr(now,30);
        OUTPUT_DIR = ['./results/','ugrid_spatial_temporal_conv', '_', timestamp,'/']
        mkdir(OUTPUT_DIR);

        frep = fopen([OUTPUT_DIR,'report.dat'],'a');
        fprintf(frep,'%12s %5s %5s %10s %10s %10s %5s %12s %12s\n', ...
                'TOL', 'Q', 'MaxD', 'INTRP', 'CHEBIMPL', 'DT', 'TN', ...
                'InRLINF', 'OutRLINF');

        % SIMULATION PARAMETERS
        VERBOSE         = false;
        GVFREQ          = 0;
        OM              = 1;
        DIM             = 2;
        VIS             = true;

        % SPATIOAL RESOLUTION
        MAX_ERROR_PER_NODE = 1e-5; % Error per box
        RES_PER_NODE       = 10;    % Resolution per Node
        ADAPTIVE           = true;
        %INTERP_TYPE        = 'spline';
        INTERP_TYPE        = 'CHEBYSHEV';
        %CHEB_IMPL          = 'CHEBFUN';
        CHEB_IMPL          = 'IAS';
        CHEB_KIND          = 2;
        ENABLE_QMSL     = false;
        ENABLE_CQMSL    = false;
        FILTER          = true;
        TESTNUM         = 31;
        
        % ensure CQMSL => QMSL
        if ENABLE_CQMSL
            ENABLE_QMSL = true;
        end

        % TEMPORAL RESOLUTION
        TINIT   = 0;
        tfinal  = 2*pi;

        %max_level_list = [4 5 6];
        max_level_list = [3 4 5];
        tn_init = 100;
        
        MAX_LEVEL = 3;
        TN        = tn_init;
        DT        = tfinal/TN;
        advect();
        %for lvl =1:size(max_level_list,2)
        %    MAX_LEVEL = max_level_list(lvl);   % Maximum depth
        %    TN        = tn_init*2^(lvl-1);
        %    DT        = tfinal/TN;
        %    % RUN
        %    advect();
        %end
    end
end
