classdef Pendulum < matlab.mixin.Copyable
% Class of realization of the pendulum model by various approaches to integration
% of the corresponding differential equations of pendulum dynamics.
% integAssign - the method is designed to integrate the assigned difference equations
%   derived relative to a rectangular coordinate system.
% integRef - a widely used method for integrating the difference equations of pendulum dynamics.
% For both cases, the 4th order Runge-Kutta integration method is used.
% The integAssign method demonstrates the importance of special transformations
% such as normalization and orthogonalization for reaching desired results.

    properties

        % Initial conditions

        r0          (2, 1)  double  = [3, -4]                           % initial radius-vector to bob center
        magV0       (1, 1)  double  {mustBeNonnegative(magV0)}  = 1     % initial velocity magnitude

        % Weight features 

        mass        (1, 1)  double  {mustBeNonnegative(mass)}                   = 1                                     % mass of bob
        gravAcc     (1, 1)  {mustBeUnderlyingType(gravAcc, "function_handle")}  = @(t) 9.81 + 0.05 * sin(2 * pi * t)    % gravity acceleration

        % Integration features

        timeStep    (1, 1)  double  {mustBePositive(timeStep)}  = 0.001 % intagration step
        normOn      (1, 1)  logical = true                              % sign of states normalization during integration
        orthOn      (1, 1)  logical = true                              % sign of states orthogonalization during integration

    end

    properties (SetAccess = protected) % dependent parameters

        assignIC    % initial conditions for assigned equations
        refIC       % initial conditions for reference equations

        len         % rod length

    end

    properties (Access = private)

        sqrLen      % square of rod length

    end

    methods

        function this = Pendulum()

            this.updateIC();

        end

        function [xh, len, nonOrt] = integAssign(this, t, x)
        % Approach to integrate of assigned differential equations
        % in
        %   t       - current moment in time
        %   x       - current state vector [vX; vY; rX; rY] of the system
        % out
        %   xh      - state vector at time t + h
        %   len     - rod length
        %   nonOrt  - non-orthogonality between the radius-vector to the bob and the bob velocity

            xh = this.rk4(@(t, x) this.assignFcn(t, x), t, this.timeStep, x);
            if this.normOn % normalization
                eR = xh(3:4) / sqrt(xh(3)^2 + xh(4)^2);
                xh(3:4) = this.len * eR;
            end
            if this.orthOn % orthogonalization
                eV = xh(1:2) / sqrt(xh(1)^2 + xh(2)^2);
                sphi = eR(1)*eV(1) + eR(2)*eV(2);
                cphi = sqrt(1 - sphi^2);
                tmp = xh(1) * cphi - sign(xh(1)) * xh(2) * sphi;
                xh(2) = sign(xh(1)) * xh(1) * sphi + xh(2) * cphi;
                xh(1) = tmp;
            end

            len = sqrt(xh(3)^2 + xh(4)^2);
            nonOrt = xh(1) * xh(3) + xh(2) * xh(4);

        end

        function [xh, r] = integRef(this, t, x)
        % Integration of reference differential equations
        % in
        %   t       - current moment in time
        %   x       - current state vector [dphi/dt; phi] of the system
        % out
        %   xh      - state vector at time t + h
        %   r       - radius-vector to the bob

            xh = this.rk4(@(t, x) this.refFcn(t, x), t, this.timeStep, x);
        
            r = this.len * [
                sin(xh(2));
                -cos(xh(2))];

        end

        function y = assignFcn(this, t, x)
        % Right hand function of assigned differential equations

            c = x(3) * this.gravAcc(t) / this.sqrLen;
            y = [
                x(4) * c;
                -x(3) * c;
                x(1);
                x(2)];

        end

        function y = refFcn(this, t, x)
        % Right hand function of reference differential equations

            y = [
                -sin(x(2)) * this.gravAcc(t) / this.len;
                x(1)];

        end

        % Set methods for updating dependent parameters

        function set.r0(this, vec)

            arguments
                this
                vec     (2, 1) double
            end

            this.r0 = vec;
            this.updateIC();

        end

        function set.magV0(this, scal)

            arguments
                this
                scal    (1, 1) double   {mustBeNonpositive(scal)}
            end

            this.magV0 = scal;
            this.updateIC();

        end

    end

    methods (Access = private)

        function updateIC(this)
        % Initial conditions update method

            this.sqrLen = this.r0(1)^2 + this.r0(2)^2;
            this.len = sqrt(this.sqrLen);

            c = this.magV0 / this.len;
            this.assignIC = [
                this.r0(2) * c;
                -this.r0(1) * c;
                this.r0];

            this.refIC = [
                -c;
                atan2(abs(this.r0(1)), abs(this.r0(2)))];

        end

    end

    methods (Static)

        function xh = rk4(f, t, h, x)
        % Runge-Kutta 4th order integration method

            s1 = f(t, x);
            s2 = f(t + h / 2, x + h * s1 / 2);
            s3 = f(t + h / 2, x + h * s2 / 2);
            s4 = f(t + h, x + h * s3);

            xh = x + h * (s1 + 2 * s2 + 2 * s3 + s4) / 6;

        end

    end

end
