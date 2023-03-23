function sgm_flat_multislice = PCD_curve_to_flat_multislice(sgm_curve, curve_radius, curve_del)
%sgm_curve / sgm_flat: sinograms (views x detector pixels x slices)
%curve_radius: radius of curvature for curved detector (mm)
%curve_del: pixel pitch for curved detector (mm)
flat_del = curve_del;
flat_del_z = flat_del;
curve_del_z = curve_del;
%Notes: assumes detectors have even number of pixels and slices

[num_views, curve_elements, curve_slices] = size(sgm_curve);
flat_elements = curve_elements;

sgm_flat = zeros(num_views, flat_elements, curve_slices);
for view = 1:num_views
    for slice_num = 1:curve_slices
        for j = 1:flat_elements
            i = atan((j - (flat_elements+1)/2)*flat_del/curve_radius)*curve_radius/curve_del+(curve_elements+1)/2;
            i_low = floor(i);
            i_high = ceil(i);
            sgm_flat(view, j, slice_num) = sgm_curve(view, i_low, slice_num)*(i_high-i)/(i_high-i_low) + sgm_curve(view, i_high, slice_num)*(i-i_low)/(i_high-i_low);
        end
    end
end

% generate multislice flat panel sinogram with same number of slices as
% original curved sinogram
sgm_flat_multislice = zeros(num_views, flat_elements, curve_slices);
for view = 1:num_views
    for k_prime = 1:curve_slices
        for j = 1:flat_elements
            i = atan((j - (flat_elements+1)/2)*flat_del/curve_radius)*curve_radius/curve_del+(curve_elements+1)/2;
            k = (curve_slices+1)/2 + (k_prime-(curve_slices+1)/2)*flat_del_z/curve_del_z*cos(curve_del/curve_radius*(i-(curve_elements+1)/2));
            k_low = floor(k);
            k_high = ceil(k);
            sgm_flat_multislice(view, j, k_prime) = sgm_flat(view, j, k_low)*(k_high-k)/(k_high-k_low) + sgm_flat(view, j, k_high)*(k-k_low)/(k_high-k_low);
        end
    end 
end