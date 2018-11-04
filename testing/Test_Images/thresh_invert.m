function im_stack_new = thresh_invert(im_stack,bit)

thresh = prctile(im_stack(:),95);

im_stack(im_stack > thresh) = thresh;

im_stack_new = 2^bit-1-im_stack;%inverts