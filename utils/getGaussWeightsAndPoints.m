function [weights, gaussPoints] = getGaussWeightsAndPoints(order)
    % Check if the order is valid
    if order < 1 || order > 14
        error('Invalid order for Gauss integration.');
    end

    % Initialize the matrices with zeros
    weights = zeros(order, 1);
    gaussPoints = zeros(order, 1);

    if order == 1
        weights(1) = 2.0;
        gaussPoints(1) = 0.0;
    elseif order == 2
        weights(1) = 1.0;
        weights(2) = 1.0;
        gaussPoints(1) = -0.577350269189626;
        gaussPoints(2) = 0.577350269189626;
    elseif order == 3
        weights(1) = 0.555555555555556;
        weights(2) = 0.888888888888889;
        weights(3) = 0.555555555555556;
        gaussPoints(1) = -0.774596669241483;
        gaussPoints(2) = 0.0;
        gaussPoints(3) = 0.774596669241483;
    elseif order == 14
        % 14-point Gauss integration rule.
        cor = 0.335180055401662;
        cen = 0.886426592797784;
        
        G14 = [-cor, -cor, -cor, -cor, +cor, +cor, +cor, +cor, -cen, cen, 0.0, 0.0, 0.0, 0.0;
               -cor, -cor, +cor, +cor, -cor, -cor, +cor, +cor, 0.0, 0.0, -cen, cen, 0.0, 0.0;
               -cor, +cor, -cor, +cor, -cor, +cor, -cor, +cor, 0.0, 0.0, 0.0, 0.0, -cen, cen];
        
        W14 = zeros(order, 1);
        W14(1:8) = cor;
        W14(9:14) = cen;

        weights = W14;
        gaussPoints = G14;
    end
end
