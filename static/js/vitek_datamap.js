$(document).ready(function () {
  console.log("loading data vitek");
  var a = "{{defaultvalues}}";
  console.log(a);
  $(".submit_data").click(function () {
    //selected objects
    console.log(a);

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
      index_value.push($(this).text());
    });
    var index_values_str = index_value.toString();
    // column value
    var column_value = [];
    $("#sortable2 li").each(function () {
      column_value.push($(this).text());
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
      $("#pivotable").html("");
      if (response["msg"]) {
        saveData(response["table"], "pivottable.csv");
        $("#pivotable").append(response["msg"])
        var data = JSON.parse(response["table_tofront"])
        console.log(typeof (data));
        var data1 = [{ "astatus": 0, "acreated_at": "2023-01-06T04:27:39.511Z", "aupdated_at": null, "adeleted_at": null, "acreated": "orgdb", "aupdated": null, "adeleted": null, "orgbatch_id": "GN_0751_01", "card_type": "AST", "card_code": "AST-GN96", "expiry_date": "2022-10-08", "instrument": "00001A0FD535 (IMB Vitek)", "proc_date": "2022-02-25", "analysis_time": "9.08 hours" }, { "astatus": 0, "acreated_at": "2023-01-06T04:27:39.426Z", "aupdated_at": null, "adeleted_at": null, "acreated": "orgdb", "aupdated": null, "adeleted": null, "orgbatch_id": "GN_0750_01", "card_type": "AST", "card_code": "AST-GN96", "expiry_date": "2022-10-08", "instrument": "00001A0FD535 (IMB Vitek)", "proc_date": "2022-02-25", "analysis_time": "11.00 hours" }, { "astatus": 0, "acreated_at": "2023-01-06T04:27:39.338Z", "aupdated_at": null, "adeleted_at": null, "acreated": "orgdb", "aupdated": null, "adeleted": null, "orgbatch_id": "GN_0749_01", "card_type": "AST", "card_code": "AST-GN96", "expiry_date": "2022-10-08", "instrument": "00001A0FD535 (IMB Vitek)", "proc_date": "2022-02-25", "analysis_time": "8.60 hours" }]
        create_pivottable(data)
      } else {
        data = response["table"];
        $("#pivotable").append(data);
      }
    })
    .fail((XMLHttpRequest, textStatus, errorThrown) => {
      console.log(XMLHttpRequest, textStatus, errorThrown);
    });
};

function create_pivottable(data) {
  $('#output').pivot(
    data, {
    cols: ['expiry_date'],
    rows: ["analysis_time"],
    aggregatorName: "intSum",
    vals: ["proc_date"],
    rendererName: "Table"
  }
  )
}