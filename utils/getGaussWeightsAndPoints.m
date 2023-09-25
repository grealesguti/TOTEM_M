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
    elseif order == 4
        % 4-point Gauss integration rule.
        weights = [0.6521451548625461; 0.6521451548625461; 0.3478548451374538; 0.3478548451374538];
        gaussPoints = [-0.1834346424956498; 0.1834346424956498; -0.5255324099163290; 0.5255324099163290];
    elseif order == 5
        % 5-point Gauss integration rule.
        weights = [0.5688888888888889; 0.4786286704993665; 0.4786286704993665; 0.2369268850561891; 0.2369268850561891];
        gaussPoints = [0.0000000000000000; -0.5384693101056831; 0.5384693101056831; -0.9061798459386640; 0.9061798459386640];
    elseif order == 14
        % 14-point Gauss integration rule.
          cor=0.758786910639328;
          cen=0.795822425754222;
          G14=[ -cor,-cor,-cor,-cor,+cor,+cor,+cor,+cor,    -cen,cen,0,0,0,0
                -cor,-cor,+cor,+cor,-cor,-cor,+cor,+cor,    0,0,-cen,cen,0,0    
                -cor,+cor,-cor,+cor,-cor,+cor,-cor,+cor,    0,0,0,0,-cen,cen
          ];
          cor=0.335180055401662;
          cen=0.886426592797784;
          W14(1:8)=cor; W14(9:14)=cen;
        weights = W14;
        gaussPoints = G14;
    end
end
