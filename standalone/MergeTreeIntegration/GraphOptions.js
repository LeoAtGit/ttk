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
    d3.select("#topN").on("input", e => {
      const topN = e.target.value;
      this.callback_on_input_change(topN);
    });

    const topN = d3.select("#topN").node().value;
    this.callback_on_input_change(topN)
  }

  callback_on_input_change(topN) {
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

    d3.selectAll("#color_legend .entry").remove();

    labels_array.forEach(d => {
      let entry = d3.select("#color_legend")
        .append("div")
        .attr("class", "entry")

      entry.append("div")
        .attr("style", `background-color: ${d[1]}; width: 10%; float: left;`)
        .text("O");  // so the div is shown...

      entry.append("div")
        .attr("style", `width: 90%; float: right;`)
        .text(this.donutGraph.labels[d[0]]);
    });
  }
}