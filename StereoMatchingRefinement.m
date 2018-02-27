function [x1,y1,z1,t1] = StereoMatchingRefinement(H1, H2 ,matches, I1, I2)

iter = 1;
maxiter = 10;
dx = 100;

u = [matches(:,1:2) ones(size(matches,1),1)];
v = [matches(:,6:7) ones(size(matches,1),1)];
H1 = H1/H1(3,3);
H2 = H2/H2(3,3);
x = [H1(2,:)';H1(3,:)';H2(2,:)';H2(3,:)'];

while norm(dx) > 1e-2 && iter < maxiter
        % visualization
    transform_matches_im1 = [matches(:, 1:2) ones(size(matches,1), 1)];
    transform_matches_im2 = [matches(:, 6:7) ones(size(matches,1), 1)];
    
    H1_i = [H1(1,:)*x(6);x(1:3)';x(4:6)'];
    H2_i = [H2(1,:)*x(12);x(7:9)';x(10:12)'];
    
    transform_matches_im1 = H1_i*transform_matches_im1';
    transform_matches_im1(1,:) = transform_matches_im1(1,:)./transform_matches_im1(3,:);
    transform_matches_im1(2,:) = transform_matches_im1(2,:)./transform_matches_im1(3,:);

    transform_matches_im2 = H2_i*transform_matches_im2';
    transform_matches_im2(1,:) = transform_matches_im2(1,:)./transform_matches_im2(3,:);
    transform_matches_im2(2,:) = transform_matches_im2(2,:)./transform_matches_im2(3,:);
    
    transform_matches = [transform_matches_im1(1:2,:)' matches(:, 3:5) transform_matches_im2(1:2,:)' matches(:, 8:10)];
    
    tform1 = projective2d(H1_i');
    ov = imref2d([size(I1) 3]);
    im1 = imwarp(I1, tform1, 'OutputView', ov);

    tform2 = projective2d(H2_i');
    im2 = imwarp(I2, tform2, 'OutputView', ov);

    VisualizeMatches(im1, im2, transform_matches, 'b', 'g');
    pause;
    
    
    % building Jacobian and residual
    u = [matches(:, 1:2) ones(size(matches,1),1)];
    v = [matches(:, 6:7) ones(size(matches,1),1)];

    J = zeros(size(matches,1), 12);
    res = zeros(size(matches,1), 1);
    for i = 1:size(matches,1)
       J(i,1:6) = 1/(u(i,:)*x(4:6))*[u(i,:) (-u(i,:)*x(1:3)/(u(i,:)*x(4:6)))*u(i,:)];
       J(i,7:12) = 1/(v(i,:)*x(10:12))*[-v(i,:) (v(i,:)*x(7:9)/(v(i,:)*x(10:12)))*v(i,:)];
       res(i) =  u(i,:)*x(1:3)/(u(i,:)*x(4:6)) - v(i,:)*x(7:9)/(v(i,:)*x(10:12));       
       % test jacobian:
%        dx = 0.01*randn(12,1);
%        xbar = x + dx;
%        
%        f = u(i,:)*xbar(1:3)/(u(i,:)*xbar(4:6)) - v(i,:)*xbar(7:9)/(v(i,:)*xbar(10:12))
%        f_hat = res(i)
%        f_hat_J = res(i) + J(i,:)*dx       
    end
        
    norm(res)
    
    % dx = (J'*J + 1000000*eye(12))\(J'*res);
    % x = x + dx;    
    
    
    % Descent direction
    descent_dir = -J'*res*1e-6;
    step_size = 0.5;
    beta = 0.1;
    alpha = 0.5;
    f = norm(res); f_p = norm(res) + 10;
    % Armijo step size:
    % loop until find the one smaller than alpha sloped linear prediction
    while((f_p - f > alpha * step_size * (J'*res)'*descent_dir))
        step_size = step_size * beta; 
        assert(step_size > 1e-16, 'step size is too small');

        x_p = x + step_size *descent_dir;
        
        for i = 1:size(matches,1)            
           res(i) =  u(i,:)*x_p(1:3)/(u(i,:)*x_p(4:6)) - v(i,:)*x_p(7:9)/(v(i,:)*x_p(10:12));       
        end  
        
        f_p = norm(res)
        alpha * step_size * (J'*res)'*descent_dir
    end
    f_p
    x = x_p;
    
end