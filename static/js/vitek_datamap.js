var savecsv = null;
$(document).ready(function () {
  console.log("loading data vitek");

  $(".submit_data").click(function () {
    var selected_data = [];
    $("input:checkbox[name=type]:checked").each(function () {
      selected_data.push($(this).val().toString());
    });
    var card_barcode = $("input[name=card_barcode]").val();
    // plot value
    var value_str = $("[data-name=data_process_value] option:selected")
      .val()
      .toString();
    console.log(value_str);
    // plot functijon
    var data_function_str = $("[data-name=data_function_type] option:selected")
      .val()
      .toString();

    // index value
    var index_value = [];
    $("#sortable3 li").each(function () {
      index_value.push($(this).data('name'));
    });
    var index_values_str = index_value.toString();
    // column value
    var column_value = [];
    $("#sortable2 li").each(function () {
      column_value.push($(this).data('name'));
    });
    var column_value_str = column_value.toString();
    // all data
    var data = {
      selected_data: selected_data,
      values: value_str,
      columns: column_value_str,
      index: index_values_str,
      card_barcode: card_barcode,
      functions: data_function_str,
    };

    sendToServer(data);
  });
});

const csrftoken = getCookie("csrftoken");
const sendToServer = (data) => {
  console.log("send to server");
  // console.log(data);

  $.ajax({
    url: "/vitekcard_list", //url,
    type: "POST",
    headers: { "X-CSRFToken": csrftoken },
    data: data,
  })
    .done((response) => {
      console.log(response)
      data = response["table_html"];
      savecsv = response["table_csv"];
      // console.log(savecsv)
      $("#pivotable").html("");
      $("#pivotable").append(data);
    })
    .fail((XMLHttpRequest, textStatus, errorThrown) => {
      console.log(XMLHttpRequest, textStatus, errorThrown);
    });
};

$("#download_pt_as_csv").click(function () {
  console.log("clicked!");
  saveData(savecsv, "pivottable.csv");
});
// function create_pivottable(data) {
//   var index_value = [];
//   $("#sortable3 li").each(function () {
//     index_value.push($(this).text());
//   });
//   var column_value = [];
//   $("#sortable2 li").each(function () {
//     column_value.push($(this).text());
//   });
//   var value_str = $("[data-name=data_process_value] option:selected")
//     .val()
//     .toString();
//   $("#output").pivot(data, {
//     cols: column_value,
//     rows: index_value,
//     aggregatorName: "intSum",
//     vals: [value_str],
//     rendererName: "Table",
//   });
// }
//
