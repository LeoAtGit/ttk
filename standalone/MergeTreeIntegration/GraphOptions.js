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
    // TopN
    d3.select("#topN").on("input", d => {
      const topN = d.target.value;
      d3.select("#label_topN").text(`${topN} TopN Categories`);

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

      d3.select("#color_legend")
        .selectAll("div")
        .data(labels_array)
        .join(
          enter => enter
            .append("div")
            .attr("style", d => `background-color: ${d[1]}`)
            .text(d => this.donutGraph.labels[d[0]]),
          update => update
            .attr("style", d => `background-color: ${d[1]}`)
            .text(d => this.donutGraph.labels[d[0]]),
          exit => exit.remove(),
        );
    });

    d3.select("#label_topN").text(`${d3.select("#topN").node().value} TopN Categories`);

    // d3.select("#color_legend").append()
  }
}