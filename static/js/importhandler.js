// custom javascript
$(document).ready(() => {
  console.log("document ready!");
});
const csrftoken = getCookie("csrftoken");
if ($("#filepath").text()) {
  $("#progressbar span:first-child").toggleClass("bg-success");
}
// Function to call the function run_task in Views.py

$(".button").on("click", function () {
  $("#preLoader").fadeIn();
  console.log("button clicked");
  var filepathlist = [];
  $("input[name=uploadedfiles_select]:checked").each(function () {
    filepathlist.push($(this).val().toString());
  });

  const datamodel = $("#datamodel").text() ? $("#datamodel").text() : "none";

  $.ajax({
    url: "/import-VITEK/",
    data: {
      type: $(this).data("type"),
      filepathlist: filepathlist,
      datamodel: datamodel.toString(),
    },
    method: "POST",
    headers: { "X-CSRFToken": csrftoken },
  })
    .done((res) => {
      console.log(res.table_name);
      const html = `
      <tr>
      <td>${res.table_name}</td>
      <td>${res.validate_result}</td>
      <td class= "text-truncate">${res.file_report}</td>
      </tr>`;
      $("#tasks").append(html);
      if (res.validate_result.includes("True")) {
        $("#Import_step3").toggleClass("visible");
        $("#progressbar span:nth-child(2)").toggleClass("bg-success");
      } else {
        $("#save_Proceed").prop("disabled", true);
        $("#save_Proceed").addClass("disabled");
        $("#next_to_confirm").toggleClass("visible");
      }

      $("#preLoader").fadeOut();
    })
    .fail((err) => {
      console.log(err);
    });
});

$(".confirmButton").on("click", function () {
  $("#preLoader").fadeIn();
  var filepathlist = [];
  $("input[name=uploadedfiles_select]:checked").each(function () {
    filepathlist.push($(this).val().toString());
  });
  const datamodel = $("#datamodel").text() ? $("#datamodel").text() : "none";
  $.ajax({
    url: "/import-VITEK/",
    method: "POST",
    data: {
      type: $(this).data("type"),
      filepathlist: filepathlist,
      datamodel: datamodel.toString(),
    },
    headers: { "X-CSRFToken": csrftoken },
  })
    .done((res) => {
      $("#preLoader").fadeOut();
      console.log(res);
      // if (res.task_status === "Form is Valid") {
      $("#Import_step4").toggleClass("visible");
      $("#progressbar span:nth-child(3)").toggleClass("bg-success");
      // }
      const html = `<p>${res.status}</p>`;
      $("#mesg_save_Proceed").append(html);
      // save_data(res);
    })
    .fail((err) => {
      console.log(err);
    });
});
