class GraphOptions {
  constructor() {
  }

  setDonutGraph(donut) {
    this.donutGraph = donut;
  }

  setStreamGraph(stream) {
    this.streamGraph = stream;
  }

  render() {
    d3.select("#color_scheme").on("input", e => {
      this.streamGraph.changeColorScheme(e.target.value);
      this.streamGraph.render("streamgraph");
      this.streamGraph.doTransform();

      this.donutGraph.changeColorScheme(e.target.value);
      this.donutGraph.render("donut");
      this.donutGraph.doTransform();

      // so the color legend will be updated as well
      d3.select("#topN").node().dispatchEvent(new Event("input"));
    });

    d3.select("#topN").on("input", e => {
      const topN = e.target.value;
      d3.select("#label_topN").text(`${topN}`);

      this.streamGraph.streamgraph_options.topN = parseInt(topN);
      this.streamGraph.render("streamgraph");
      this.streamGraph.doTransform();

      this.donutGraph.donut_options.topN = parseInt(topN);
      this.donutGraph.render("donut");
      this.donutGraph.doTransform();

      // render the color legend
      let labels = this.donutGraph.tree.getLabels();

      // transform the labels to an array that is sorted by color_order
      let labels_array = [];
      for (let key in labels) {
        labels_array.push([key, labels[key]])
      }
      labels_array = labels_array.sort((a, b) =>
        this.donutGraph.tree.color_order.indexOf(parseInt(a[0])) - this.donutGraph.tree.color_order.indexOf(parseInt(b[0])));

      d3.selectAll("#color_legend .entry").remove();

      labels_array.forEach(d => {
        let entry = d3.select("#color_legend")
          .append("div")
          .attr("class", "entry row")

        entry.append("div")
          .attr("style", `background-color: ${d[1]};`)
          .attr("class", "col-1 rounded-circle");

        entry.append("div")
          .attr("class", "col-10")
          .text(this.donutGraph.labels[d[0]]);
      });
    });

    d3.select("#maxwidthRoot").on("input", e => {
      const maxwidthRoot = e.target.value;
      d3.select("#label_maxwidthRoot").text(`${maxwidthRoot}`);

      this.streamGraph.streamgraph_options.maxwidth_root = parseInt(maxwidthRoot);
      this.streamGraph.resetLayoutingCoords();
      this.streamGraph.render("streamgraph");
      this.streamGraph.doTransform();
    });

    d3.select("#absolute_width").on("change", e => {
      const absolute_width = e.target.checked;

      this.streamGraph.streamgraph_options.use_relative_sizes = !absolute_width;
      this.streamGraph.streamgraph_options.maxwidth_root = parseInt(d3.select("#maxwidthRoot").node().value);
      this.streamGraph.resetLayoutingCoords();
      this.streamGraph.render("streamgraph");
      this.streamGraph.doTransform();
    });

    d3.select("#contour_count").on("input", e => {
      kdeRenderer.update_nContours(parseInt(e.target.value));
      d3.select("#label_contour_count").text(`${e.target.value}`);
    });

    // so the values are displayed when the side is loaded
    d3.select("#maxwidthRoot").node().dispatchEvent(new Event("input"));
    d3.select("#absolute_width").node().dispatchEvent(new Event("change"));
    d3.select("#topN").node().dispatchEvent(new Event("input"));
    d3.select("#color_scheme").node().dispatchEvent(new Event("input"));
    d3.select("#contour_count").node().dispatchEvent(new Event("input"));
  }
}