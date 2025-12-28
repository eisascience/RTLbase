# scoped seeds enable snapshot-friendly artifacts

    Code
      seeded_df
    Output
             x         y
    1 0.2655087 0.8983897
    2 0.3721239 0.9446753
    3 0.5728534 0.6607978
    4 0.9082078 0.6291140
    5 0.2016819 0.0617863

    Code
      ggplot2::ggplot_build(plot_obj)$data[[1]][c("x", "y")]
    Output
             x         y
    1 0.2655087 0.8983897
    2 0.3721239 0.9446753
    3 0.5728534 0.6607978
    4 0.9082078 0.6291140
    5 0.2016819 0.0617863

