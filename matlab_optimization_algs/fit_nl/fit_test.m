function fit_test(which)
% This function tests the fit_nl and fit_nl_ex routines. Several fits in 1D-3D
% with built-in functions or custom defined functions are performed.
% 
% Syntax
%   function fit_test()
%   function fit_test(which)
%
% Input parameter
%   which   which test 1-7 (or nothing for a menu)

if nargin == 0

    % option which not given, print a menu
    which = 1;
    while which ~= 0
        % the menu
        fprintf('Which test would You like to do?\n');
        fprintf(' (1) Calling built-in function in fit_func()\n');
        fprintf(' (2) 1D fit of noisy data with a Lorentzian function\n');
        fprintf(' (3) 2D fit of noisy data with a custom defined function\n');
        fprintf(' (4) 3D fit of noisy data with a Gaussian function\n');
        fprintf(' (5) 1D fit of noisy data with 3 Lorentzian peaks\n');
        fprintf(' (6) 2D fit of noisy data with 4 Gaussian peaks and complex side conditions\n');
        fprintf(' (7) Runtime comparison of lsqcurvefit(), fmincon() and fminsearch()\n');
        fprintf(' (Enter) Exit\n');
        which = input('Choose [Enter]:', 's');
        if isempty(which)
            % which is empty when enter was pressed
            which = 0;
        else
            which = int32(str2double(which));
            if which > 0 && which < 8
                % call the right test
                fit_test(which);
            else
                %  all other cases (not only Enter) we also exit
                which = 0;
            end
        end
    end
else

    % according to which execute a test
    fprintf('\n');
    switch which
        case 1
            % a test: Calling built-in functions in fit_func() ------------

            disp('Calling built-in function in fit_func()');

            % x values run from - 1 to 5
            xdata.x = -1:0.1:5;

            % function parameter for function 1, a gaussian with amplitude
            % 10, center at 0,5, FWHM 2 and background 0.1
            x = [0.1, 10, 0.5, 2];
            % compute it using a built-in function
            y1 = fit_func('1d_gaussian', x, xdata);

            % using the same x-values we have the parameters for function
            % 2, a sine with amplitude 4, cycle duration 3, start at -1
            % and background of 3
            x = [3., 4., 3, -1];
            y2 = fit_func(@custom_func_1dsine, x, xdata);

            % plot both in one graph
            figure;
            hold on;
            plot(xdata.x, y1, 'b');
            plot(xdata.x, y2, 'r');
            title('function data created by a call to fit\_func()');

        case 2
            % tests 2-4 use function fit_nl()
            % a test: 1D fitting with a Lorentzian ------------------------

            disp('1D fit of noisy data with a Lorentzian function');

            % prepare image and initial parameters
            xdata.x = 0.5:99.5; % 100 x-values between 0 and 100

            % the 1D Lorentzian model has 4 parameters: background,
            % amplitude, center and FWHM
            x_true = [0.5, 20., 35., 10.];

            % left and right borders
            lb   = [0., 1., 0., 0.1];
            ub   = [1., 100., 100., 20.];

            % calculate the function corresponding to x_true
            y_true = fit_func('1d_lorentzian', x_true, xdata);

            % add gaussian noise to y_true to obtain y_meas
            y_meas = random('Normal', y_true, 3); % standard deviation = 3

            % set reasonable start values for fitting parameters
            x0 = [0.1, 10., 50., 5.];

            % fit y_meas with start values x0 with a call to fit_nl
            % parameter fixed (5th) is not set
            x_fit = fit_nl('1d_lorentzian', x0, xdata, y_meas, [], lb, ub);

            % calculate y_fit corresponding to x_fit
            y_fit = fit_func('1d_lorentzian', x_fit, xdata);

            % plot everything in a figure
            figure;
            hold on;
            plot(xdata.x, y_true, 'b');
            plot(xdata.x, y_meas, 'k');
            plot(xdata.x, y_fit, 'r');
            ylim([-10 25]);
            title('test 1 - 1D fitting of a Lorentzian peak');
            text(2, -8, 'blue - original curve, black - noisy curve, red - fitted curve');

            % display parameters overview
            disp('original parameters');
            disp(x_true);
            disp('fitted parameters');
            disp(x_fit);

        case 3
            % test 2D fitting of a sine with a custom curve ---------------

            % Comment: Here we use a custom function and give two fields x
            % and y in xdata so that the function knows on which grid in
            % the R^2 (the 2D world) the function should be computed. x and
            % y have the same size as the output and can be conveniently by
            % created by a call to ndgrid(). However, every other way,
            % where xdata is a string, a cell array or whatsever, which
            % makes your function know whose x and y values (of the
            % independent variables) to take is also okay.

            disp('2D fit of noisy data with a custom defined function');

            % prepare xdata / computational grid
            [x, y]  = ndgrid(0.5:19.5, 0.5:19.5);
            xdata.x = x;
            xdata.y = y;
            % our custom function has 4 parameters (see custom_func_2dsine
            % for a function expression)

            % the true parameter
            x_true = [10., 0.1, 0.1, 0.2];

            % calculate the 2D function corresponding to x_true
            y_true = fit_func(@custom_func_2dsine, x_true, xdata);

            % add gaussian noise to y_true to obtain y_meas
            y_meas = random('Normal', y_true, 6); % standard deviation = 6

            % set reasonable start values for fitting parameters
            x0 = [8., 0.12, 0.12, 0];

            % fit y_meas with start values x0 with a call to fit_nl
            % parameter fixed, lb and ub are not set
            x_fit = fit_nl(@custom_func_2dsine, x0, xdata, y_meas);

            % calculate y_fit corresponding to x_fit
            y_fit = fit_func(@custom_func_2dsine, x_fit, xdata);

            % plot everything in a figure
            clims = [-100 220];
            subplot2(2, 2, 1);
            imagesc(y_true, clims);
            title('original image');
            subplot2(2, 2, 2);
            imagesc(y_meas, clims);
            title('noisy image');
            subplot2(2, 2, 3);
            imagesc(y_fit, clims);
            title('fitted image');
            subplot2(2, 2, 4);
            imagesc(y_true - y_fit, clims);
            title('difference original - fitted image');

            % display parameters overview
            disp('original parameters');
            disp(x_true);
            disp('fitted parameters');
            disp(x_fit);

        case 4
            % test 3D fitting of a gaussian -------------------------------

            disp('3D fit of noisy data with a Gaussian function');

            % prepare xdata / computational grid
            [x, y, z]  = ndgrid(0.5:19.5, 0.5:19.5, 0.5:19.5);
            % in all directions 20 points ranging vom 0.5 to 19.5
            xdata.x = x;
            xdata.y = y;
            xdata.z = z;

            % the 3D Gaussian function has 7 parameters (see fit_func.m)

            % the true parameters
            x_true = [0.1, 50., 8., 11., 9., 2., 1.5, 3.];

            % calculate the 2D function corresponding to x_true
            y_true = fit_func('3d_gaussian', x_true, xdata);

            % add gaussian noise to y_true to obtain y_meas
            y_meas = random('Normal', y_true, 6); % standard deviation = 6

            % set reasonable start values for fitting parameters
            x0 = [0.5, 47., 7., 12., 8.5, 3., 1.8, 2.7];

            % the first parameter (background) shall be kept fixed
            fixed = [1, 0, 0, 0, 0, 0, 0, 0];

            % fit y_meas with start values x0 with a call to fit_nl
            % parameter lb and ub are not set
            [x_fit resnorm] = fit_nl('3d_gaussian', x0, xdata, y_meas, fixed);

            % calculate y_fit corresponding to x_fit
            y_fit = fit_func('3d_gaussian', x_fit, xdata);

            % we do not plot (would need to be a 3D plot)
            fprintf('3D fit completed with fixed parameter 1, residual norm is %.2f\n', resnorm);

            % display parameters overview
            disp('original parameters');
            disp(x_true);
            disp('initial guess for parameters');
            disp(x0);
            disp('fitted parameters');
            disp(x_fit);

        case 5
            % test 5 and 6 use function fit_nl_ex()
            % test: 1D fit of 3 Lorentzian peaks --------------------------

            disp('1D fit of noisy data with 3 Lorentzian peaks');

            % prepare xdata / computational grid
            xdata = [];
            xdata.x  = 0.5:99.5;

            % Three 1D Lorentzian functions have 1+3*(2*1+1)=10 parameter
            % (see fitfunc.m)

            % the true parameters
            x_true = [0.1, 20., 20., 8., 10., 50., 5., 15., 70., 6.5];

            % calculate the 2D function corresponding to x_true
            y_true = fit_func('1d_lorentzian', x_true, xdata);

            % add gaussian noise to y_true to obtain y_meas
            y_meas = random('Normal', y_true, 3); % standard deviation = 3

            % set reasonable start values for fitting parameters
            x0 = [0.3, 18., 23., 7., 12., 55., 3., 17., 75., 4.5];

            % create options
            options.opt = optimset('Display', 'iter');

            % fit y_meas with start values x0 with a call to fit_nl_ex
            x_fit = fit_nl_ex('1d_lorentzian', x0, xdata, y_meas, options);

            % since we did not do any special, the fit did exactly the same
            % as with fit_nl, but this demonstrates only, that fmincon and
            % lsqcurvefit relies on the same principles

            % calculate y_fit corresponding to x_fit
            y_fit = fit_func('1d_lorentzian', x_fit, xdata);

            % plot everything in a figure
            figure;
            hold on;
            plot(xdata.x, y_true, 'b');
            plot(xdata.x, y_meas, 'k');
            plot(xdata.x, y_fit, 'r');
            ylim([-10 25]);
            title('test 5 - 1D fitting of 3 Lorentzian peaks');
            text(2, -4, 'blue - original curve, black - noisy curve, red - fitted curve');
            text(2, -6, 'watch the center peak which is sometimes not fitted correctly');
            text(2, -8, 'because of noise and a bad starting point');

            % display parameters overview
            disp('original parameters');
            disp(x_true);
            disp('fitted parameters');
            disp(x_fit);

        case 6
            % test: 2D of 4 gaussian peaks --------------------------------

            disp('2D fit of noisy data with 4 Gaussian peaks and complex side conditions');
            % now our showcase :)))

            % prepare xdata / computational grid
            xdata = [];
            [x, y] = ndgrid(0.5:99.5, 0.5:99.5);
            xdata.x = x;
            xdata.y = y;

            % 4 2D Gaussian functions have 1+4*(2*2+1)=21 parameter
            % thats a lot (see fitfunc.m)

            % the true parameters
            x_true = 0.1;   % background
            x_true = [x_true, 20., 30., 30., 8., 10.];  % first peak
            x_true = [x_true, 22., 70., 60., 6., 9.];   % second peak
            x_true = [x_true, 26., 40., 70., 11., 7.];  % third peak
            x_true = [x_true, 20., 80., 90., 8., 8.];   % fourd peak

            % calculate the 2D function corresponding to x_true
            y_true = fit_func('2d_gaussian', x_true, xdata);

            % add gaussian noise to y_true to obtain y_meas
            y_meas = random('Poisson', y_true); % standard deviation = sqrt(mean)

            % set reasonable start values for fitting parameters
            x0 = 0.3; % background
            x0 = [x0, 15., 25., 35., 8.5, 8.5];     % first peak
            x0 = [x0, 20., 65., 65., 8.5, 8.5];     % second peak
            x0 = [x0, 30., 45., 65., 8.5, 8.5];     % third peak
            x0 = [x0, 15., 75., 85., 8.5, 8.5];     % fourd peak

            % calculate the 2D function corresponding to x0
            y0 = fit_func('2d_gaussian', x0, xdata);

            % create options
            options.lb = zeros(size(x0)); % all parameters should be at least positive
            options.ub = zeros(size(x0)) + 100; % and at most 100
            options.likelihood  = 'Poisson';
            options.penalty_fun = @penalty_fun;

            % fit y_meas with start values x0 with a call to fit_nl_ex
            fprintf('Fitting in progress, please wait!');
            x_fit = fit_nl_ex('2d_gaussian', x0, xdata, y_meas, options);

            % compute image corresponding to x_fit
            y_fit = fit_func('2d_gaussian', x_fit, xdata);

            % plot everything in a figure
            clims = [0 30];
            subplot2(2, 2, 1);
            imagesc(y_true, clims);
            colormap(hot);
            title('original image');
            subplot2(2, 2, 2);
            imagesc(y_meas, clims);
            title('noisy image');
            subplot2(2, 2, 3);
            imagesc(y0, clims);
            title('initial guess image');
            subplot2(2, 2, 4);
            imagesc(y_fit, clims);
            title('fitted image');

            % display parameters overview
            fprintf('parameters\n');
            fprintf('true\t initial\t fitted\n');
            for ki = 1:length(x_true)
                fprintf('%.1f\t %.1f\t %.1f\n', x_true(ki), x0(ki), x_fit(ki));
            end

        case 7
            % test: Runtime comparison lsqcurvefit/fmincon/fminsearch -----

            disp('Measure running time for a simple example: 1D Gaussian peak');
            % prepare a simple example - a 1D gaussian peak
            
            % All default options for maxIter, TolX, etc. should be the same
            % for lsqcurvefit(), fmincon(), fminsearch()
            % A good sign for this is that the reached residual values were
            % the same in my case

            % prepare image and initial parameters
            xdata.x = 0.5:99.5; % 100 x-values between 0 and 100

            % the 1D Gaussian model has 4 parameters: background,
            % amplitude, center and FWHM
            x_true = [0.5, 20., 35., 10.];

            % calculate the function corresponding to x_true
            y_true = fit_func('1d_gaussian', x_true, xdata);

            % to have more runs to add we create several noisy copies of
            % our data and fit them all
            N = 100;
            fprintf('Performing %d fitting starts per functions.\n', N);

            % add gaussian noise to y_true to obtain y_meas
            y_true2 = repmat(y_true, N, 1);
            y_meas2 = random('Normal', y_true2, 3); % standard deviation = 3

            % set reasonable start values for fitting parameters
            x0 = [0.1, 15., 40., 8.];

            % call to lsqcurvefit via fit_nl
            fprintf('Fitting via lsqcurvefit\n');
            tic
            r1 = 0;
            for ki = 1 : N
                [x r] = fit_nl('1d_gaussian', x0, xdata, y_meas2(ki, :));
                r1 = r1 + r;
            end
            r1 = r1 / N;
            toc
            fprintf('Mean residual value (resnorm of lsqcurvefit) %.2f\n', r1);

            % call to fmincon via fit_nl_ex
            fprintf('Fitting via fmincon\n');
            tic
            r2 = 0;
            for ki = 1 : N
                [x r] = fit_nl_ex('1d_gaussian', x0, xdata, y_meas2(ki, :));
                r2 = r2 + r;
            end
            r2 = r2 / N;
            toc
            fprintf('Mean residual value (fval of fmincon) %.2f\n', r2);

            % call to fminsearch via fit_nl_ex
            options = [];
            options.use_fmincon = 'no';
            fprintf('Fitting via fminsearch\n');
            tic
            r3 = 0;
            for ki = 1 : N
                [x r] = fit_nl_ex('1d_gaussian', x0, xdata, y_meas2(ki, :), options);
                r3 = r3 + r;
            end
            r3 = r3 / N;
            toc
            fprintf('Mean residual value (fval of fminsearch) %.2f\n', r3);

        otherwise
            error('Variable which out of range!');
    end
    fprintf('\n');

end

% local functions ---------------------------------------------------------

% custom function for test 1
    function y = custom_func_1dsine(x, xdata)
        % has 4 parameters and expects xdata.x to hold the x values
        % form: p1 + p2 * sin(((x-p4)/p3)*2pi)
        y = zeros(size(xdata.x)) + x(1);
        y = y + x(2) * sin((xdata.x - x(4)) / x(3) * 2 * pi);
    end

% custom function for test 3
    function y = custom_func_2dsine(x, xdata)
        % is 2D (expects xdata.x and xdata.y)
        % form: p1 * sin(p2 * x) * cos(p3 * y) + p4
        y = x(1) * sin(x(2) * xdata.x) * cos(x(3) * xdata.y) + x(4);
    end

% penalty function for test 6
    function y = penalty_fun(x)
        % simply adds a quadratic penalty if amplitude of the 4 gaussian peaks
        % is far away from 22, e.g. this would help against overfitting with
        % darker (relative to true brightness) peaks if the number of peaks
        % would be unknown, a weighting factor should be introduced
        % additionally
        y = (x(2)-22)^2+(x(7)-22)^2+(x(12)-22)^2+(x(17)-22)^2;
    end

end