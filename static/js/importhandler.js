// custom javascript
$(document).ready(() => {
  console.log("document ready!");
});
var output_files = $("#output");

$("#id_file_field").change(function () {
  output_files.empty();
  var files = $(this).prop("files");
  var htmls = ``;
  jQuery.each(files, function (i, val) {
    htmls += `<p>${i + 1}-${val.name}</p>`;
  });
  output_files.append(htmls);
});
const csrftoken = getCookie("csrftoken");
if ($("#filepath").text()) {
  $("#progressbar span:first-child").toggleClass("bg-success");
}
// Function to call the function run_task in Views.py

$(".button").on("click", function () {
  $("#preLoader").fadeIn();

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
      if (res.status === "Delete") {
        $("#delete_cancel_task a:first-child").addClass("blink");
      } else {
        var validateResult = JSON.parse(res.validate_result.replace(/'/g, '"'));
        var validateReport = JSON.parse(
          res.file_report
            .replace(/'/g, '"')
            .replace(/"{/g, "{")
            .replace(/}"/g, "}")
        );
        console.log(validateReport);
        // JSON.parse;
        var f_list = Object.keys(validateResult);
        // ------------
        for (let i = 0; i < f_list.length; i++) {
          var error_num = 0;
          var warning_num = 0;
          var ew_description = { Error: [], Warning: [] };
          validateReport[f_list[i]].forEach((el) => {
            if (el["Error"] > 0) {
              error_num++;
              ew_description["Error"].push(el["Error"].toString());
            }
            if (el["Warning"] > 0) {
              warning_num++;
              ew_description["Warning"].push(el["Warning"].toString());
            }
          });

          const tr = `
          <tr>
          <td>${f_list[i]}</td>
          <td>${validateResult[f_list[i]]}</td>
          <td>${error_num.toString()}</td>
      <td>${warning_num.toString()}</td>
      <td>${JSON.stringify(ew_description)}</td>
      </tr>`;
          $("#tasks").append(tr);

          if (error_num === 0) {
            $("#progressbar span:nth-child(2)").toggleClass("bg-success");
          } else {
            $(".confirmButton").prop("disabled", true);
          }
        }
        // --------------
        $("#tasksreport").append(`
    <div id="upload_report">${res.file_report}</div>`);

        $("#Import_step3").addClass("visible");

        $("#preLoader").fadeOut();
      }
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
