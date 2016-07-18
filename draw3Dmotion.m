function [] = draw3DmotionV2(t, x, y, z, phi, theta, psi)
    
    x_max = max(abs(x));
    y_max = max(abs(y));
    z_max = max(abs(z));
    
    side_max = max(x_max, y_max);
    side_max = max(side_max, 155);
    
    l = side_max/3;
    
    rotated_corner1 = ones(3, length(t));
    rotated_corner2 = ones(3, length(t));
    rotated_corner3 = ones(3, length(t));
    rotated_corner4 = ones(3, length(t));
    
    
    %the unrotated corners are at (l, -l, 0), (l, l, 0),
    %(-l, l, 0), (-l, -l, 0)
    for i = 1:length(t)
        rotation_matrix = ...
[cos(psi(i))*cos(phi(i))-sin(psi(i))*cos(theta(i))*sin(phi(i)), -sin(psi(i))*cos(phi(i))-cos(psi(i))*cos(theta(i))*sin(phi(i)), sin(theta(i))*sin(phi(i));
 cos(psi(i))*sin(phi(i))+sin(psi(i))*cos(theta(i))*cos(phi(i)), -sin(psi(i))*sin(phi(i))+cos(psi(i))*cos(theta(i))*cos(phi(i)), -sin(theta(i))*cos(phi(i));
 sin(psi(i))*sin(theta(i)), cos(psi(i))*sin(theta(i)), cos(theta(i))];
       
        rotated_corner1(:,i) = rotation_matrix*[l/sqrt(2); l/sqrt(2); 0];
                            
        rotated_corner2(:,i) = rotation_matrix*[-l/sqrt(2); l/sqrt(2); 0];
                        
        rotated_corner3(:,i) = rotation_matrix*[-l/sqrt(2); -l/sqrt(2); 0];
                        
        rotated_corner4(:,i) = rotation_matrix*[l/sqrt(2); -l/sqrt(2); 0];
    end
    
    i = 1;
    h = plot3(x(1:i), y(1:i), z(1:i), ...
        x(i)+rotated_corner1(1, i), y(i)+rotated_corner1(2, i), z(i)+rotated_corner1(3, i), '*', ...
        x(i)+rotated_corner2(1, i), y(i)+rotated_corner2(2, i), z(i)+rotated_corner2(3, i), '*', ...
        x(i)+rotated_corner3(1, i), y(i)+rotated_corner3(2, i), z(i)+rotated_corner3(3, i), '*', ...
        x(i)+rotated_corner4(1, i), y(i)+rotated_corner4(2, i), z(i)+rotated_corner4(3, i), '*');
    
    for i = 1:length(x)
        set(h, {'XData'}, {x(i); x(i)+rotated_corner1(1, i);  x(i)+rotated_corner2(1, i); x(i)+rotated_corner3(1, i); x(i)+rotated_corner4(1, i)});
        set(h, {'YData'}, {y(i); y(i)+rotated_corner1(2, i);  y(i)+rotated_corner2(2, i); y(i)+rotated_corner3(2, i); y(i)+rotated_corner4(2, i)});
        set(h, {'ZData'}, {z(i); z(i)+rotated_corner1(3, i);  z(i)+rotated_corner2(3, i); z(i)+rotated_corner3(3, i); z(i)+rotated_corner4(3, i)});
        xlabel('x'); ylabel('y'); zlabel('z');
        title(['Time: ', num2str(t(i))]);
        %h = plot3(xp(1:i), yp(1:i), zp(1:i), rotated_corner1(1, i), rotated_corner1(2, i), rotated_corner1(3, i), '*');
        grid on;
        xlim([-1.5*side_max 1.5*side_max]);
        ylim([-1.5*side_max 1.5*side_max]);
        zlim([-5*side_max 5*side_max]);
        drawnow;
    end    
end