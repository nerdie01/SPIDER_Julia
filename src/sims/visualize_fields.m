function visualize_fields(f1, f2, titlef1, titlef2)
    figure(1);
    clf;

    subplot(1,2,1);
    imagesc(real(ifft2(f1)));
    axis square;
    colorbar();
    title(titlef1);

    subplot(1,2,2);
    imagesc(real(ifft2(f2)));
    axis square;
    colorbar();
    title(titlef2);

    drawnow;
end
