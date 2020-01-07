# Interactive-Visualization-with-R

ukb dataset is the example dataset.

## IBD Segments Plot

```{bash}
king -b ukb.bed --ibdseg --prefix ukb
```
Output files:
ukb.seg,ukb.segments.gz,ukballsegs.txt


```{bash}
Run IBDseg.R in the same directory.
```

Reference: http://people.virginia.edu/~wc9c/KING/manual.html#IBDSEG

## ROH segments plot

```{bash}
king -b ukb.bed --roh --prefix ukb
```
Output files:
ukballsegs.txt,ukb.roh,ukb.rohseg.gz


```{bash}
Run ROH.R in the same directory.
```

Reference: http://people.virginia.edu/~wc9c/KING/KINGvisualization.html



