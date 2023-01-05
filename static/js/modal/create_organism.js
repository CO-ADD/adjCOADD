$(document).ready(function () {
  let my_modal = $("#createOrganismModal");

  $("#createOrganism").click(function () {
    my_modal.load("/createOrg/", function (response) {
      if (response == "Permission Not Granted") {
        alert("permission!");
      } else {
        my_modal.modal("show");
      }
    });
  });
});
