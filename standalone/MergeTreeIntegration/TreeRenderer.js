class TreeRenderer {
  constructor(treeContainer){

    this.vtkDataSet = null;
    this.treeContainer = treeContainer;

    this.treeContainer.append("<svg></svg>");
    this.svg = d3.select("svg")
        .attr("height", 500)
        .attr("width", 500);

    this.root = this.svg.append("g");
    this.nodelayer = this.root.append("g");
  }

  setVtkDataSet(dataset){
    this.vtkDataSet = dataset;
  }

  render(){
    if (!this.vtkDataSet)
	  console.error("trying to render without dataset");

    this.nodelayer.empty();
    const coords = this.vtkDataSet.points.coordinates.data;
    const branchid = this.vtkDataSet.pointData.BranchId.data;
    const npoints = coords.length / 3;
    const scale = 150;

    for (let i = 0; i < npoints; i++)
    {
      this.nodelayer.append("circle")
          .attr("cx", coords[i*3] * scale)
          .attr("cy", coords[i*3 + 1] * scale)
          .attr("r", 3)
          .attr("fill", branchid[i] === 0 ? "red" : "green");
    }
  }
}

