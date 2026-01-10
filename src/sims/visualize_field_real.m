function visualize_field_real(f1, titlef1)
    figure(1);
    clf;

    imagesc(f1);
    axis square;
    colorbar();
    title(titlef1);

    drawnow;
end