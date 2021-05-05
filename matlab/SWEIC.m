function q0 = SWEIC(x, y)
    hu0 = 0;
    hv0 = 0;
%     if (x > 0 && y > 0)
%         h0 = 1;
%     else
%         h0 = 0.5;
%     end


%     A = 2;
%     x0 = 0;
%     y0 = 0;
%     sigma_x = 0.5;
%     sigma_y = 0.5;
%     h0 = A*exp(-(((x - x0)^2)/(2*sigma_x^2) + ((y - y0)^2)/(2*sigma_y^2)));

    h0 = 0.5;
    q0 = [h0; hu0; hv0];
    
end