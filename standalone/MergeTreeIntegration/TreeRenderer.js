class TreeRenderer {
  constructor(treeContainer){
    this.treeContainer = treeContainer;

    this.scale = 200;
    this.padding = 50;
    this.width = 500;
    this.height = 500;

    this.treeContainer.append("<svg></svg>");
    this.svg = d3.select("svg")
        .attr("height", this.width)
        .attr("width", this.height)
        .attr("style", "border:1px solid black");

    this.root = this.svg.append("g");
    this.nodelayer = this.root.append("g");

    this.svg.call(d3.zoom()
        // .extent([[0, 0], [500, 500]])
        // .scaleExtent([1, 8])
        .on("zoom", ({transform}) => this.root.attr("transform", transform))
    );

    this.vtkDataSet = null;
    this.points = [];
  }

  setVtkDataSet(dataset){
    this.vtkDataSet = dataset;

    // point information
    const _coords = this.vtkDataSet.points.coordinates.data;
    let coords = Array.from(_coords);  // convert from Float64Array to normal Array, so we can use `coords.shift()`
    const branchid = this.vtkDataSet.pointData.BranchId.data;
    const kde = this.vtkDataSet.pointData.KDE.data;
    const kde_i0 = this.vtkDataSet.pointData.KDE_ByType_I0.data;
    const kde_i1 = this.vtkDataSet.pointData.KDE_ByType_I1.data;

    // change the y coordinates, because ttk and svg use a different coordinate system
    let _y = [];
    for (let i = 1; i < coords.length; i+=3) {
      _y.push(coords[i]);
    }
    const _y_max = Math.max(..._y);
    const _y_min = Math.min(..._y);

    while(coords.length > 0) {
      this.points.push(new Point(
          coords.shift() * this.scale,
          ((Math.abs(_y_max - _y_min)) - (coords.shift() - _y_min)) * this.scale,
          coords.shift() * this.scale)
      );
    }

    for(let i = 0; i < this.points.length; i++) {
      this.points[i].setBranchId(branchid[i]);
      this.points[i].setKDE(kde[i]);
      this.points[i].setKDE_I0(kde_i0[i]);
      this.points[i].setKDE_I1(kde_i1[i]);
    }

    // center the tree with viewBox
    const x_max = Math.max(...this.points.map(p => p.x));
    const x_min = Math.min(...this.points.map(p => p.x));
    const y_max = Math.max(...this.points.map(p => p.y));
    const y_min = Math.min(...this.points.map(p => p.y));

    this.svg
        .attr("viewBox", [
            x_min - this.padding,
            y_min - this.padding,
            2 * this.padding + (x_max - x_min),
            2 * this.padding + (y_max - y_min)]
        );

    // cell information
    this.connectivityArray = this.vtkDataSet.cells.connectivityArray.data;
    const offsetsArray = this.vtkDataSet.cells.offsetsArray.data;

    // sanity check that we only work with lines, and that the data for
    // the lines is 'nicely' arranged in the connectivityArray.
    //
    // If there is connectivity information for other geometric primitives,
    // this has to be implemented separately.
    for (let i = 0; i < offsetsArray.length; i++) {
      if (BigInt(i*2) !== offsetsArray[i]) {
        console.error("there were non-line segments in the data");
        console.error("alternatively, the data was not 'nicely' arranged (i.e. always i*2)");
        return;
      }
    }
  }

  render() {
    if (!this.vtkDataSet) {
      console.error("trying to render without dataset");
      return;
    }

    this.nodelayer.empty();

    // draw the lines of the graph (in an unoptimized way)
    const line = d3.line()
        .x(d => d.x)
        .y(d => d.y);

    for (let i = 0; i < this.connectivityArray.length; i+=2) {
      this.nodelayer.append("path")
          .attr("fill", "none")
          .attr("stroke", "steelblue")
          .attr("stroke-width", 1.5)
          .attr("stroke-miterlimit", 1)
          .attr("d", line([this.points[this.connectivityArray[i]], this.points[this.connectivityArray[i+1]]]));
    }

    for (let i = 0; i < this.points.length; i++) {
      this.nodelayer.append("circle")
          .attr("cx", this.points[i].x)
          .attr("cy", this.points[i].y)
          .attr("r", 3)
          .attr("fill", this.points[i].branchId === 0 ? "red" : "green");
    }
  }
}

class Point {
  constructor(x, y, z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  setBranchId(id) {
    this.branchId = id;
  }

  setKDE(kde) {
    this.kde = kde;
  }

  setKDE_I0(kde) {
    this.kde_i0 = kde;
  }

  setKDE_I1(kde) {
    this.kde_i1 = kde;
  }
}
