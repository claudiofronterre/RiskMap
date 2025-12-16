# Plot Spatial Scores for a Specific Model and Metric

This function visualizes spatial scores for a specified model and
metric. It combines test set data, handles duplicate locations by
averaging scores, and creates a customizable map using ggplot2.

## Usage

``` r
plot_score(object, which_score, which_model, ...)
```

## Arguments

- object:

  A list containing test sets and model scores. The structure should
  include \`object\$test_set\` (list of sf objects) and
  \`object\$model\[\[which_model\]\]\$score\[\[which_score\]\]\`.

- which_score:

  A string specifying the score to visualize. Must match a score
  computed in the model.

- which_model:

  A string specifying the model whose scores to visualize.

- ...:

  Additional arguments to customize ggplot, such as
  \`scale_color_gradient\` or \`scale_color_manual\`.

## Value

A ggplot object visualizing the spatial distribution of the specified
score.
