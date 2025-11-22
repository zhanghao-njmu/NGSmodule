/**
 * Confirm Dialog Component - Reusable confirmation dialogs
 * 统一的确认对话框组件
 */
import { Modal } from 'antd'
import { ExclamationCircleFilled, DeleteOutlined, WarningOutlined } from '@ant-design/icons'

export interface ConfirmDialogOptions {
  title?: string
  content?: string | React.ReactNode
  onConfirm?: () => void | Promise<void>
  onCancel?: () => void
  okText?: string
  cancelText?: string
  okButtonProps?: any
  type?: 'info' | 'success' | 'warning' | 'error'
}

/**
 * Generic confirm dialog
 */
export const confirm = ({
  title = '确认操作',
  content,
  onConfirm,
  onCancel,
  okText = '确定',
  cancelText = '取消',
  okButtonProps,
  type = 'warning',
}: ConfirmDialogOptions) => {
  Modal.confirm({
    title,
    content,
    okText,
    cancelText,
    okButtonProps,
    icon: <ExclamationCircleFilled style={{ color: '#f59e0b' }} />,
    onOk: onConfirm,
    onCancel,
    centered: true,
  })
}

/**
 * Delete confirmation dialog
 */
export const confirmDelete = (
  itemName: string,
  onConfirm: () => void | Promise<void>,
  options?: Partial<ConfirmDialogOptions>
) => {
  Modal.confirm({
    title: '确认删除',
    content: (
      <div>
        <p>您确定要删除 <strong>{itemName}</strong> 吗？</p>
        <p style={{ color: '#ef4444', fontSize: 13 }}>此操作不可撤销，所有相关数据将被永久删除。</p>
      </div>
    ),
    icon: <DeleteOutlined style={{ color: '#ef4444' }} />,
    okText: '删除',
    cancelText: '取消',
    okButtonProps: { danger: true },
    onOk: onConfirm,
    centered: true,
    ...options,
  })
}

/**
 * Batch delete confirmation dialog
 */
export const confirmBatchDelete = (
  count: number,
  onConfirm: () => void | Promise<void>,
  options?: Partial<ConfirmDialogOptions>
) => {
  Modal.confirm({
    title: '批量删除确认',
    content: (
      <div>
        <p>您确定要删除选中的 <strong>{count}</strong> 项吗？</p>
        <p style={{ color: '#ef4444', fontSize: 13 }}>此操作不可撤销，所有相关数据将被永久删除。</p>
      </div>
    ),
    icon: <DeleteOutlined style={{ color: '#ef4444' }} />,
    okText: `删除 ${count} 项`,
    cancelText: '取消',
    okButtonProps: { danger: true },
    onOk: onConfirm,
    centered: true,
    ...options,
  })
}

/**
 * Dangerous action confirmation dialog
 */
export const confirmDangerousAction = (
  action: string,
  description: string,
  onConfirm: () => void | Promise<void>,
  options?: Partial<ConfirmDialogOptions>
) => {
  Modal.confirm({
    title: `确认${action}`,
    content: (
      <div>
        <p>{description}</p>
        <p style={{ color: '#ef4444', fontSize: 13, marginTop: 8 }}>
          ⚠️ 警告：此操作具有危险性，请谨慎操作！
        </p>
      </div>
    ),
    icon: <WarningOutlined style={{ color: '#ef4444' }} />,
    okText: action,
    cancelText: '取消',
    okButtonProps: { danger: true },
    onOk: onConfirm,
    centered: true,
    ...options,
  })
}

/**
 * Leave page confirmation dialog
 */
export const confirmLeave = (
  onConfirm: () => void | Promise<void>,
  message = '您有未保存的更改，确定要离开吗？'
) => {
  Modal.confirm({
    title: '未保存的更改',
    content: message,
    icon: <ExclamationCircleFilled style={{ color: '#f59e0b' }} />,
    okText: '离开',
    cancelText: '继续编辑',
    okButtonProps: { danger: true },
    onOk: onConfirm,
    centered: true,
  })
}

/**
 * Task execution confirmation dialog
 */
export const confirmTaskExecution = (
  taskName: string,
  estimatedTime: string,
  onConfirm: () => void | Promise<void>
) => {
  Modal.confirm({
    title: '确认执行任务',
    content: (
      <div>
        <p>任务名称：<strong>{taskName}</strong></p>
        <p>预计耗时：<strong>{estimatedTime}</strong></p>
        <p style={{ marginTop: 12, color: '#6b7280' }}>
          任务开始后将在后台运行，您可以在任务列表中查看进度。
        </p>
      </div>
    ),
    icon: <ExclamationCircleFilled style={{ color: '#2563eb' }} />,
    okText: '开始执行',
    cancelText: '取消',
    onOk: onConfirm,
    centered: true,
  })
}

/**
 * Pipeline execution confirmation dialog
 */
export const confirmPipelineExecution = (
  pipelineName: string,
  sampleCount: number,
  onConfirm: () => void | Promise<void>
) => {
  Modal.confirm({
    title: '确认运行管道',
    content: (
      <div>
        <p>管道名称：<strong>{pipelineName}</strong></p>
        <p>样本数量：<strong>{sampleCount}</strong></p>
        <p style={{ marginTop: 12, color: '#6b7280' }}>
          管道将对所有样本执行完整的分析流程，这可能需要较长时间。
        </p>
        <p style={{ color: '#f59e0b', fontSize: 13 }}>
          请确保已正确配置所有参数和输入文件。
        </p>
      </div>
    ),
    icon: <ExclamationCircleFilled style={{ color: '#2563eb' }} />,
    okText: '开始运行',
    cancelText: '检查参数',
    onOk: onConfirm,
    centered: true,
  })
}

export default {
  confirm,
  confirmDelete,
  confirmBatchDelete,
  confirmDangerousAction,
  confirmLeave,
  confirmTaskExecution,
  confirmPipelineExecution,
}
