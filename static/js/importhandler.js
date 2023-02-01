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
      var file_list = $.grep(
        res.validatefile_name.split(","),
        (n) => n == 0 || n
      );

      var validateResult = JSON.parse(res.validate_result.replace(/'/g, '"'));
      // var validateReport = JSON.parse(
      //   res.file_report.replace(/'/g, '"').replace(/\\/g, "")
      // );
      // JSON.parse;
      for (let i = 0; i < file_list.length; i++) {
        const tr = `
      <tr>
      <td>${file_list[i]}</td>
      <td>${validateResult[file_list[i]]}</td>
      </tr>`;
        $("#tasks").append(tr);
      }
      $("#tasksreport").append(`
        <div id="upload_report">${res.file_report}</div>`);

      $("#Import_step3").addClass("visible");
      if (res.validate_result.includes("True")) {
        $("#progressbar span:nth-child(2)").toggleClass("bg-success");
      } else {
        $(".confirmButton").prop("disabled", true);
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
      $("#Import_step4").addClass("visible");
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
