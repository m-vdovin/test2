classdef PendulumTest < matlab.unittest.TestCase
% Class for testing the Pendulum class integration methods

    properties

        tspan       = [0, 100]                                  % integration bounds
        tolInteg    = 1e-4                                      % integration tolerance

        natOdeTol   = odeset('RelTol', 1e-6, 'AbsTol', 1e-7)    % native ODE solver tolerance settings

    end

    methods (Test)

        function integrationTest(this)

            obj     = Pendulum;
            time    = this.tspan(1):obj.timeStep:this.tspan(end);
            sol     = [];
            ref     = [];

            sol(:, 1) = time;
            x0 = obj.assignIC;
            sol(1, 2:5) = x0;
            sol(1, 6) = obj.len;
            sol(1, 7) = x0(1) * x0(3) + x0(2) * x0(4);
            assignX = x0;

            ref(:, 1) = time;
            x0 = obj.refIC;
            ref(1, 2:3) = x0;
            ref(1, 4) = obj.len * sin(x0(2));
            ref(1, 5) = -obj.len * cos(x0(2));
            refX = x0;

            for idx = 2:size(time, 2)
                [assignX, sol(idx, 6), sol(idx, 7)] = obj.integAssign(time(idx), assignX);
                sol(idx, 2:5) = assignX;

                [refX, ref(idx, 4:5)] = obj.integRef(time(idx), refX);
                ref(idx, 2:3) = refX;
            end

            accuracyInteg = max(abs(sol(:, 4:5) - ref(:, 4:5)), [], "all");
            this.verifyLessThan(accuracyInteg, this.tolInteg);

            nat = deval(ode45(@(t, x) obj.refFcn(t, x), this.tspan, x0, this.natOdeTol), time)';
            accuracyInteg = max(abs(ref(:, 2:3) - nat(:, 1:2)), [], "all");
            this.verifyLessThan(accuracyInteg, this.tolInteg);

        end

    end

end
