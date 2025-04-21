function fsa_quantity = flux_surface_average(quantity, jacobian, theta, phi, volume)
    if nargin < 5
        volume = false; % Default value for volume if not provided
    end

    % Initialize the integrals
    denom = 0;
    numerator = 0;

    % Number of points in theta and phi
    nTheta = length(theta);
    nPhi = length(phi);

    % Compute the denominator
    for i = 1:nPhi
        for j = 1:nTheta
            if i == 1
                % Left edge
                if j == 1
                    denom = denom + (phi(1,i+1) - phi(1,i)) * (theta(1,j+1) - theta(1,j)) * jacobian(j, i);
                elseif j == nTheta
                    % Right edge
                    denom = denom + (phi(1,i+1) - phi(1,i)) * (theta(1,j) - theta(1,j-1)) * jacobian(j, i);
                else
                    % Middle points
                    denom = denom + (phi(1,i+1) - phi(1,i)) * (theta(1,j+1) - theta(1,j-1))/2 * jacobian(j, i);
                end
            elseif i == nPhi
                % Last phi edge
                if j == 1
                    denom = denom + (phi(1,i) - phi(1,i-1)) * (theta(1,j+1) - theta(1,j)) * jacobian(j, i);
                elseif j == nTheta
                    denom = denom + (phi(1,i) - phi(1,i-1)) * (theta(1,j) - theta(1,j-1)) * jacobian(j, i);
                else
                    denom = denom + (phi(1,i) - phi(1,i-1)) * (theta(1,j+1) - theta(1,j-1))/2 * jacobian(j, i);
                end
            else
                % Middle points
                if j == 1
                    denom = denom + (phi(1,i+1) - phi(1,i-1))/2 * (theta(1,j+1) - theta(1,j)) * jacobian(j, i);
                elseif j == nTheta
                    denom = denom + (phi(1,i+1) - phi(1,i-1))/2 * (theta(1,j) - theta(1,j-1)) * jacobian(j, i);
                else
                    denom = denom + (phi(1,i+1) - phi(1,i-1))/2 * (theta(1,j+1) - theta(1,j-1))/2 * jacobian(j, i);
                end
            end
        end
    end

    % Compute the numerator
    for i = 1:nPhi
        for j = 1:nTheta
            if i == 1
                % Left edge
                if j == 1
                    numerator = numerator + (phi(1,i+1) - phi(1,i)) * (theta(1,j+1) - theta(1,j)) * quantity(j, i) * jacobian(j, i);
                elseif j == nTheta
                    % Right edge
                    numerator = numerator + (phi(1,i+1) - phi(1,i)) * (theta(1,j) - theta(1,j-1)) * quantity(j, i) * jacobian(j, i);
                else
                    % Middle points
                    numerator = numerator + (phi(1,i+1) - phi(1,i)) * (theta(1,j+1) - theta(1,j-1))/2 * quantity(j, i) * jacobian(j, i);
                end
            elseif i == nPhi
                % Last phi edge
                if j == 1
                    numerator = numerator + (phi(1,i) - phi(1,i-1)) * (theta(1,j+1) - theta(1,j)) * quantity(j, i) * jacobian(j, i);
                elseif j == nTheta
                    numerator = numerator + (phi(1,i) - phi(1,i-1)) * (theta(1,j) - theta(1,j-1)) * quantity(j, i) * jacobian(j, i);
                else
                    numerator = numerator + (phi(1,i) - phi(1,i-1)) * (theta(1,j+1) - theta(1,j-1))/2 * quantity(j, i) * jacobian(j, i);
                end
            else
                % Middle points
                if j == 1
                    numerator = numerator + (phi(1,i+1) - phi(1,i-1))/2 * (theta(1,j+1) - theta(1,j))  * quantity(j, i) * jacobian(j, i);
                elseif j == nTheta
                    numerator = numerator + (phi(1,i+1) - phi(1,i-1))/2 * (theta(1,j) - theta(1,j-1))  * quantity(j, i) * jacobian(j, i);
                else
                    numerator = numerator + (phi(1,i+1) - phi(1,i-1))/2 * (theta(1,j+1) - theta(1,j-1))/2  * quantity(j, i) * jacobian(j, i);
                end
            end
        end
    end

    % Calculate flux surface average quantity
    fsa_quantity = numerator / denom;

    if volume
        % Return the denominator if volume is true
        fsa_quantity = denom;
    end
end