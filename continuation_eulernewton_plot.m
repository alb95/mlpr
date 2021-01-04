function [x, it, did_turnaround] = continuation_eulernewton_plot(target_alpha, v, R, tol, maxit)

% modeled on alg. 6.1.10 in [Georg, Allgower] by Federico Poloni and Alberto Bucci

deltatilde = 20;

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

n = length(v);
I = eye(n);
it = 0;
did_turnaround = false;

% get a starting point
alpha = 0.6;
[x, itn] = newton(alpha, v, R, tol, maxit - it);
fprintf('Initial fixed-alpha Newton required %d iterations:\n', itn);
it = it + itn;
y = [x; alpha];

previous_direction = [zeros(n,1);1];
h = 0.01;

clf;
hold on;
plot(alpha, x(1), 'xb');

while true
    it = it + 1; %we count 1 for each predictor _or_ corrector iteration
    y_old = y;
    x = y(1:end-1); alpha = y(end);
    JH = [alpha*R*(kron(x, I) + kron(I, x)) - I, R*kron(x, x)-v];
    [Q, ~] = qr(JH');
    t = Q(:, end); % already normalized with norm(t)==1
    % makes sure we are moving "in the right direction"
    % Georg-Allgower use determinants instead, but they are not easily
    % available with Matlab's QR.
    if previous_direction' * t < 0
        t = -t;
    end
    % predictor
    y = y_old + h*t;
    x = y(1:end-1); alpha = y(end);
    % adjust h so that we never pass target_alpha (little point in doing
    % that)
%     if alpha > target_alpha
%        h = (target_alpha - y_old(end)) / t(end);
%        y = y_old + h*t;
%        x = y(1:end-1); alpha = y(end);
%     end
    % plot([y_old(end) alpha], [y_old(1) x(1)], '-r');
    quiver(y_old(end), y_old(1), alpha-y_old(end), x(1)-y_old(1), 'r', 'MaxHeadSize', 0.01);
    H = alpha*R*kron(x,x) + (1-alpha)*v - x;
    % corrector
    y_before_corrector = y;
    k = 0;
    repeat_predictor = false;
    while norm(H, 1) > tol
        k = k + 1;
        it = it + 1;
        JH = [alpha*R*(kron(x, I) + kron(I, x)) - I, R*kron(x, x)-v];
        [Q, r] = qr(JH');
        Q = Q(:,1:n);
        r = r(1:n, :);
        newton_step = - Q*(r'\H); %Moore-Penrose inverse is Q*inv(r');
        if k == 1
            delta = norm(newton_step,1) / h^2;
            f = sqrt(delta / deltatilde);
            if f > 2
                repeat_predictor = true;
                break;
            end
        end
        y = y + newton_step; 
        x = y(1:end-1); alpha = y(end);
        H = alpha*R*kron(x,x) + (1-alpha)*v - x;
        if k + it >= maxit
            return
        end
    end
    if repeat_predictor
        fprintf('Newton step was too long: repeating predictor with smaller step size\n');
        h = h/2;
        continue;
    end
    fprintf('Corrector computed solution with alpha=%g and residual %g in %d iterations:\n', alpha, norm(H,1), k);
    plot(alpha, x(1), 'xb');
    if y(end) < y_old(end)
        did_turnaround = true;
    end
    if alpha > target_alpha
        break
    end
    % update step size
    f = max(f, 1/2); f = min(f, 2); %clamps f to [1/2, 2]
    h = h / f;
    h = min(h, 0.1); % limits max step size
    % set up iteration
    previous_direction = y - y_old;
end
% y_old and y contain two solutions that bracket target_alpha;
% we interpolate between them
t = y - y_old;
h = (target_alpha - y_old(end)) / t(end);
y = y_old + h*t;
x = y(1:end-1); alpha = y(end);
assert(abs(alpha - target_alpha) < tol);
[x, itn] = newton(target_alpha, v, R, tol, maxit - it, x);
fprintf('Final fixed-alpha Newton required %d iterations:\n', itn);
it = it + itn;
