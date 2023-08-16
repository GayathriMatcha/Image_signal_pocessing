function H = finding_h(points, corresponding_points, epsilon, fraction, n_samples, n_iters)
    % number of points
    n_points = length(points);
    % square now so that we do not have to take square root again and again
    epsilon = epsilon^2;
    % d as defined in the assignment PDF
    d = ceil(fraction*n_points);
    
    % Get the set of 4 points, their correspondence points, and the
    % remaining points and correspondence point

    % Get all the indices
    n_index = 1:n_points;
    % randomly select n_samples indices from all the indices
    select_indices = randperm(n_points, n_samples);
    % Get the points and corresponding points for those indices
    pts = points(select_indices, :);
    corres_pts = corresponding_points(select_indices, :);
    % Get the remaining indices using set differentiation
    % For example, if all_indices is [1, 2, 3, 4, 5] and 
    % sample_indices is [1,3], setdiff(all_indices, sample_indices);
    % gives [2, 5, 5]
    remaining_indices = setdiff(n_index, select_indices);
    % return the remaining points too to for the consensus set
    rem_pts = points(remaining_indices, :);
    rem_corres_pts = corresponding_points(remaining_indices, :);
    %% 
    % Build the A matrix
    %A = A_matrix(pts, corres_pts);
    A = zeros(2*n_samples, 9);
    for i = 1:n_samples
        x = pts(i, 1);
        y = pts(i, 2);
        x_prime = corres_pts(i, 1);
        y_prime = corres_pts(i, 2);
        % Fill the rows of A with the appropriate values
        A(2*i-1, :) = [-x -y -1  0  0  0 x*x_prime y*x_prime x_prime];
        A(2*i, :)   = [ 0  0  0 -x -y -1 x*y_prime y*y_prime y_prime];
    end
%% 
    % Compute H from A
    % Initialize the value of H
    %H = compute_homography(A);
    
    [~, ~, v] = svd(A);    % Perform SVD on the A matrix
    h = v(:, end);         % The last column of v is the solution for h
    H = reshape(h, 3, 3)'; % reshape h appropriately to get the (3,3) H matrix
%%     
    % Estimate corresponding points using the remaining points (those not
    % used in the calculation of H) and the H estimated
    %rem_pts_prime = corresponding_points_2D_multiple(H, rem_pts);
    % rem_points is of the shape (n,2) with rem_points(i,:) = [x(i), y(i)]
    % points is not in homogenous coordinates
    % add a row of 1s to bring points to homogenous coordinates
    % homogenous_points(i, :) = [x(i), y(i), 1]
    homogenous_points = [rem_pts, ones(length(rem_pts), 1)];
    homogenous_points_prime = homogenous_points*H';
    % Multiply with the H matrix and perform element wise division
    % corresponding_points (i, :) = homogenous_points_prime(i,
    % 1:2)/homogenous_points_prime(i, 3)
    rem_pts_prime = homogenous_points_prime(:, 1:2)./homogenous_points_prime(:, 3);
    %% 


    % Take the squared difference and sum to get the squared errors
    squared_difference = (rem_pts_prime - rem_corres_pts).^2;
    % sum_sq_diff(i) = (x'(i) - x'''(i))^2 + (y'(i) - y'''(i))^2 
    % where x' and x''' are the notation as used in lecture 14
    sum_sq_diff = sum(squared_difference, 2);
    % Consensus set is a logical array with 1 where the sum squared
    % difference is less than epsilon^2
    consensus_set = sum_sq_diff<epsilon;
    % Get the number of elements in the consensus set
    % Initiate the n_consensus set variable
    n_consensus = sum(consensus_set);
    if n_consensus>=d
        % If the condition is met, get all the inlier and outlier points
        % using the consensus set and recompute H using all the inliers
        inlier_pts = points(consensus_set, :);
        inlier_corres_pts = corresponding_points(consensus_set, :);
        %A = A_matrix(inlier_pts, inlier_corres_pts);
        A = zeros(2*length(inlier_pts), 9);
        for i = 1:length(inlier_pts)
            x = inlier_pts(i, 1);
            y = inlier_pts(i, 2);
            x_prime = inlier_corres_pts(i, 1);
            y_prime = inlier_corres_pts(i, 2);
            % Fill the rows of A with the appropriate values
            A(2*i-1, :) = [-x -y -1  0  0  0 x*x_prime y*x_prime x_prime];
            A(2*i, :)   = [ 0  0  0 -x -y -1 x*y_prime y*y_prime y_prime];
        end
        H = compute_homography(A);
        % Escape from the function
        return;
    end
%% 
    % Go into this for loop only if the first calculated H is not
    % satisfactory
    for nth_iter = 1:n_iters
        % Same set of steps as that in the first computation of H
        [pts, corres_pts, rem_pts, rem_corres_pts] = n_rows_sampler(points, corresponding_points, n_samples);
        %A = A_matrix(pts, corres_pts);
        A = zeros(2*length(pts), 9);
        for i = 1:length(pts)
            x = pts(i, 1);
            y = pts(i, 2);
            x_prime = corres_pts(i, 1);
            y_prime = corres_pts(i, 2);
            % Fill the rows of A with the appropriate values
            A(2*i-1, :) = [-x -y -1  0  0  0 x*x_prime y*x_prime x_prime];
            A(2*i, :)   = [ 0  0  0 -x -y -1 x*y_prime y*y_prime y_prime];
        end
        H_ = compute_homography(A);

        rem_pts_prime = corresponding_points_2D_multiple(H_, rem_pts);
        squared_difference = (rem_pts_prime - rem_corres_pts).^2;
        sum_sq_diff = sum(squared_difference, 2);
        consensus_set_ = sum_sq_diff<epsilon;
        n_consensus_ = sum(consensus_set_);
        
        if n_consensus_>=n_consensus
            % If you get a better consensus set, update the consensus set
            % as well as the value of H
            n_consensus = n_consensus_;
            consensus_set = consensus_set_;
            H = H_;
        end
        
        if n_consensus>=d
            % Whenever we have a consensus set as large as we require,
            % recompute H using all the inliers
            inlier_pts = points(consensus_set, :);
            inlier_corres_pts = corresponding_points(consensus_set, :);
            %A = A_matrix(inlier_pts, inlier_corres_pts);
            A = zeros(2*length(inlier_pts), 9);
            for i = 1:length(inlier_pts)
                x = inlier_pts(i, 1);
                y = inlier_pts(i, 2);
                x_prime = inlier_corres_pts(i, 1);
                y_prime = inlier_corres_pts(i, 2);
                % Fill the rows of A with the appropriate values
                A(2*i-1, :) = [-x -y -1  0  0  0 x*x_prime y*x_prime x_prime];
                A(2*i, :)   = [ 0  0  0 -x -y -1 x*y_prime y*y_prime y_prime];
            end
            H = compute_homography(A);
            % And break from the for loop
            break;
        end
    end
    
    % If you cannot meet the consensus set, take the case where you got the
    % maximum size of the consensus set and recompute H using all the
    % inliers
    if nth_iter == n_iters
        inlier_pts = points(consensus_set, :);
        inlier_corres_pts = corresponding_points(consensus_set, :);
        %A = A_matrix(inlier_pts, inlier_corres_pts);
        A = zeros(2*length(inlier_pts), 9);
        for i = 1:length(inlier_pts)
            x = inlier_pts(i, 1);
            y = inlier_pts(i, 2);
            x_prime = inlier_corres_pts(i, 1);
            y_prime = inlier_corres_pts(i, 2);
            % Fill the rows of A with the appropriate values
            A(2*i-1, :) = [-x -y -1  0  0  0 x*x_prime y*x_prime x_prime];
            A(2*i, :)   = [ 0  0  0 -x -y -1 x*y_prime y*y_prime y_prime];
        end
        H = compute_homography(A);
    end
end