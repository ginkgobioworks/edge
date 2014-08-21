window.SequenceViewer = function ($, options) {
    var dom_id = 'sequence-viewer';
    if ('dom_id' in options) { dom_id = options.dom_id; }
    var selector = '#'+dom_id;
    $(selector).html('<table></table>');

    function clear() {
        $('table', selector).empty();
    }

    function setSequenceWithAnnotations(sequence, annotations, start_bp, rowlen) {
        function annotationToggle(ev) {
            var el = $(ev.target);
            if (el.hasClass('sequence-annotation-selected')) {
                $('tr.sequence-annotation td div').show();
                $('tr.sequence-annotation').show();
                $('tr.sequence-annotation td div.sequence-annotation-selected')
                  .removeClass('sequence-annotation-selected');
            }
            else {
                $('tr.sequence-annotation td div').hide();
                $('tr.sequence-annotation').hide();

                var id = el.data('sequence-annotation-id');
                $('tr.sequence-annotation td div.sequence-annotation-'+id)
                  .addClass('sequence-annotation-selected')
                  .show()
                  .parents('tr')
                  .show();
            }
        }

        rowlen = typeof rowlen !== 'undefined' ? rowlen : 80;
        clear();
        var rsplit = new RegExp('(.{1,'+rowlen+'})', 'g');
        var rows = sequence.match(rsplit);

        var ai = 0;
        _.map(annotations, function(a) { a.__sequence_annotation_id = ai; ai++; });

        var row_start = start_bp;
        _.map(rows, function(row) {
            var row_end = row_start+row.length-1;
            var row_annotations = _.select(annotations, function(a) {
                return (row_start <= a.base_first && a.base_first <= row_end) ||
                       (row_start <= a.base_last && a.base_last <= row_end);
            });

            var start_annotations = [];
            var end_annotations = [];

            _.map(row_annotations, function(a) {
                if (row_start <= a.base_first && a.base_first <= row_end) {
                    var i = a.base_first-row_start;
                    if (start_annotations[i] === undefined) { start_annotations[i] = []; }
                    start_annotations[i].push(a);
                }
                if (row_start <= a.base_last && a.base_last <= row_end) {
                    var i = a.base_last-row_start;
                    if (end_annotations[i] === undefined) { end_annotations[i] = []; }
                    end_annotations[i].push(a);
                }
            });

            var tr_mid = $('<tr class="sequence-bp"></tr>');
            var tr_top = [];
            var tr_top_html = '<tr class="sequence-annotation sequence-annotation-start"></tr>';
            var tr_bot = [];
            var tr_bot_html = '<tr class="sequence-annotation sequence-annotation-end"></tr>';

            var i = 0;
            _.map(row.split(''), function(bp) {
                var td = $('<td></td>');
                td.append(bp);
                tr_mid.append(td);

                if (start_annotations.length > 0) {
                    if (start_annotations[i]) {
                        _.map(start_annotations[i], function(a) {
                            var tr = $(tr_top_html);
                            if (i > 0) { tr.append('<td colspan="'+i+'"></td>'); }
                            var td = $('<td></td>');
                            var c = '&gt;';
                            if (a.strand < 0) { c = '&lt;'; }
                            var a = $('<div></div>')
                                      .addClass('sequence-annotation-'+a.__sequence_annotation_id)
                                      .data('sequence-annotation-id', a.__sequence_annotation_id)
                                      .append('['+c+a.display_name);
                            td.append(a);
                            tr.append(td);
                            if (i < row.length-1) { tr.append('<td colspan="'+(row.length-i-1)+'"></td>'); }
                            tr_top.push(tr);
                        });
                    }
                }

                if (end_annotations.length > 0) {
                    if (end_annotations[i]) {
                        _.map(end_annotations[i], function(a) {
                            var tr = $(tr_bot_html);
                            if (i > 0) { tr.append('<td colspan="'+i+'"></td>'); }
                            var td = $('<td></td>');
                            var a = $('<div></div>')
                                      .addClass('sequence-annotation-'+a.__sequence_annotation_id)
                                      .data('sequence-annotation-id', a.__sequence_annotation_id)
                                      .append(']'+a.display_name);
                            td.append(a);
                            tr.append(td);
                            if (i < row.length-1) { tr.append('<td colspan="'+(row.length-i-1)+'"></td>'); }
                            tr_bot.push(tr);
                        });
                    }
                }

                i += 1;
            });

            if (tr_top.length > 0) { tr_top[0].addClass('sequence-annotation-start-first'); }
            _.map(tr_top, function(tr) { $('table', selector).append(tr); });
            $('table', selector).append(tr_mid);
            if (tr_bot.length > 0) { tr_bot[tr_bot.length-1].addClass('sequence-annotation-end-last'); }
            _.map(tr_bot, function(tr) { $('table', selector).append(tr); });

            row_start += row.length;
        });

        $('tr.sequence-annotation td div').on('click', annotationToggle);
    }

    function setSequence(sequence, start_bp, rowlen, column) {
        rowlen = typeof rowlen !== 'undefined' ? rowlen : 80;
        column = typeof column !== 'undefined' ? column : 10;
        clear();

        var rsplit = new RegExp('(.{1,'+rowlen+'})', 'g');
        var csplit = new RegExp('(.{1,'+column+'})', 'g');
        var rows = sequence.match(rsplit);
        var tr = $('<tr></tr>');

        var s = start_bp;
        var td_s = $('<td class="sequence-pos"></td>');
        var td_d = $('<td class="sequence-data"></td>');
        var td_l = $('<td class="sequence-pos"></td>');
        for (var i=0; i<rows.length; i++) {
            var l = s+rows[i].length-1;
            td_s.append(s+'/'+(s-start_bp+1)+'<br/>');
            var columns = rows[i].match(csplit);
            for (var j=0; j<columns.length; j++) {
                td_d.append('<span>'+columns[j]+'</span>');
            }
            td_d.append('<br/>');
            td_l.append(l+'/'+(l-start_bp+1)+'<br/>');
            s = l+1;
        }
        tr.append(td_s);
        tr.append(td_d);
        tr.append(td_l);
        $('table', selector).append(tr);
    }

    return {
        clear: clear,
        setSequence: setSequence,
        setSequenceWithAnnotations: setSequenceWithAnnotations
    }
}
